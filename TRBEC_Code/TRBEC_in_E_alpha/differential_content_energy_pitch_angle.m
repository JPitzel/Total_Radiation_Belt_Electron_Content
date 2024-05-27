%given half orbit and filenames to import data from, outputs the differential content
%as a function of energy and pitch angle and the cooresponding Lstars
function [diff_content, energy, Lstar] = differential_content_energy_pitch_angle(epoch_select,satellite,orbtimes,flux_directory,DC_directory,filename_flux)

%constants
const.c = 299792458*100;%cm/s, speed of light
const.E_0 = 0.511;%MeV, electron rest energy
const.R_E = 6371*10^5;%cm, Earth radius


%import MagEIS data from given file list
flux = []; epoch = []; Lstar = []; B_eq = []; B_sat = [];
for flux_file_idx = 1:length(filename_flux)
    [data_flux_i,~] = spdfcdfread(strcat(flux_directory,filename_flux(flux_file_idx)),'Variables',{'Epoch','FEDU_CORR','FEDU_Alpha','FEDU_Energy','L_star','B_Calc','B_Eq'});
    epoch_i = data_flux_i{1};
    flux_i = data_flux_i{2}*1000;%convert 1/keV to 1/MeV
    flux_i = flux_i(:,1:end-4,:);
    alpha_sat = data_flux_i{3};
    energy = data_flux_i{4}/1000;%convert keV to MeV
    energy = energy(1:end-4);
    Lstar_i = data_flux_i{5};
    B_sat_i = data_flux_i{6}*10^-5;%convert nT to G
    B_eq_i = data_flux_i{7}*10^-5;%convert nT to G
    
    epoch = cat(1,epoch,epoch_i);
    flux = cat(3,flux,flux_i);
    Lstar = cat(1,Lstar,Lstar_i);
    B_sat = cat(1,B_sat,B_sat_i);
    B_eq = cat(1,B_eq,B_eq_i);
end

%select out data corresponding to single half orbit
[flux,epoch,Lstar,B_sat,B_eq] = half_orbit_select(flux,epoch,Lstar,B_sat,B_eq,orbtimes,epoch_select);
[alpha_length,energy_length,Lstar_length] = size(flux);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%select bounds and resolution of the interpolation
energy_min = 0.05;
energy_max = 1;
Lstar_min = 3.5;
Lstar_max = 5.0;
alpha_interp_length = 100;
energy_interp_length = round((energy_max-energy_min)/0.05 + 1);
Lstar_interp_length = round((Lstar_max-Lstar_min)/0.1 + 1);%sample L^* at every 0.1

%enable when generating differential content plots, disable for linear spacing to calculate the TRBEC
energy_logspace_toggle = 0;
plot_toggle_interp = 0;%creates plots of data and interpolation for various Lstars, just helps check interpolation is okay
plot_toggle_diff = 0;%creates plots of differential content (intermediate step of TRBEC after integrating over alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%format into a grid format (eg. alpha goes from a vector to a 3D grid) to
%be compatible with the flux data (pointwise operations are easier than for loops)
alpha_sat = repmat(alpha_sat,[1,energy_length,Lstar_length]);
energy = repmat(energy,[1,alpha_length,Lstar_length]); energy = permute(energy,[2,1,3]);
Lstar = repmat(Lstar,[1,alpha_length,energy_length]); Lstar = permute(Lstar,[2,3,1]);
B_sat = repmat(B_sat,[1,alpha_length,energy_length]); B_sat = permute(B_sat,[2,3,1]);
B_eq = repmat(B_eq,[1,alpha_length,energy_length]); B_eq = permute(B_eq,[2,3,1]);

%convert alpha_sat to alpha_eq (the formula only returns values between 0
%and 90 so the second line to to adjust to the pitch angles back to the original)
%the 90 deg pitch angle (at the satellite) splits into 2 cases (one <90 and one >90, but is taken into account below (ctrl. f: ***)
alpha_eq = asind(sqrt(((B_eq).*(sind(alpha_sat).^2))./(B_sat)));
alpha_eq(7:end,:,:) = 180 - alpha_eq(7:end,:,:);%because sind() returns a number between 0 deg and 90 deg, convert all pitch angles known to be above 90 to what they should be;
%add in flux at the loss cone, first calc alpha_eq_LC and then append to grid formatted data
h = 100*10^5;%cm because R_E is also in cm
alpha_eq_LC = asind((1 + ((1.5).*(h./const.R_E)))./(Lstar(1,1,:).^(3/2).*(4-((3./Lstar(1,1,:)).*(1+(h./const.R_E)))).^(0.25))); %equation from JF's work
alpha_eq_LC = repmat(alpha_eq_LC,[1,energy_length,1]);

%in the case that the alpha_eq_LC is greater than the lowest pitch angle energy set alpha = 0 so it doesn't interfer with the interpolation
alpha_eq_LC(alpha_eq_LC>alpha_eq(1,1,:)) = 0;

alpha_eq_LC_gt_90 = asind((1 + ((1.5).*(h./const.R_E)))./(Lstar(end,1,:).^(3/2).*(4-((3./Lstar(end,1,:)).*(1+(h./const.R_E)))).^(0.25))); %equation from JF's work
alpha_eq_LC_gt_90 = repmat(alpha_eq_LC_gt_90,[1,energy_length,1]);

%in the case that the alpha_eq_LC_gt_90 is less than the highest pitch angle energy set alpha = 180 so it doesn't interfer with the interpolation
alpha_eq_LC_gt_90(alpha_eq_LC<alpha_eq(end,1,:)) = 0;


%add in both loss cone below and above 90 deg to the data
%*** (the index 6 is specifically adding in the other case when converting between local and equatorial pitch angles)
flux = cat(1,flux(1,:,:),flux,flux(6,:,:),flux(end,:,:));
alpha_eq = cat(1,alpha_eq_LC,alpha_eq,90*ones(1,energy_length,Lstar_length),alpha_eq_LC_gt_90);
energy = cat(1,energy(1,:,:),energy,energy(6,:,:),energy(1,:,:));
Lstar = cat(1,Lstar(1,:,:),Lstar,Lstar(6,:,:),Lstar(1,:,:));


%'folds' symmetric pitch angles
alpha_eq(alpha_eq>90) = 180 - alpha_eq(alpha_eq>90);


%interpolation points
Lstar_interp = linspace(Lstar_min,Lstar_max,Lstar_interp_length);
if energy_logspace_toggle == 0
    energy_interp = linspace(energy_min,energy_max,energy_interp_length)';
end
if energy_logspace_toggle == 1
    energy_interp = logspace(log10(energy_min),log10(energy_max),energy_interp_length)';%use logspace when generating differntial content plots!!!
end

%get weighting function to interpolate later, if file of W vs alpha_eq doesn't exist, then create one
%NOTE: takes a couple seconds to generate so having a precomputed file drastically speeds it up!
try
    load 'weighting_function_vs_alpha_eq.mat' W_for_W_interp alpha_eq_for_W_interp;
catch
    [W_for_W_interp,alpha_eq_for_W_interp] = weighting_function();
    save 'weighting_function_vs_alpha_eq' W_for_W_interp alpha_eq_for_W_interp;
end

%interpolate the flux
diff_content = zeros(Lstar_interp_length,energy_interp_length);
for i = 1:Lstar_interp_length
    %select out the closest Lstar data to the interpolated data
    Lstar_slice = Lstar_interp(i);
    [~,nearest_Lstar_index] = min(abs(Lstar(1,1,:)-Lstar_slice));
    
    %define meshes for data and interpolated points
    alpha_eq_interp_at_Lstar_slice = linspace(asind((1 + ((1.5).*(h./const.R_E)))./(Lstar_slice.^(3/2).*(4-((3./Lstar_slice).*(1+(h./const.R_E)))).^(0.25))),90,alpha_interp_length);
    %alpha_eq_interp_at_Lstar_slice = linspace(51.2172,90,alpha_interp_length);
    
    [energy_interp_at_Lstar_slice,alpha_eq_interp_at_Lstar_slice] = meshgrid(energy_interp,alpha_eq_interp_at_Lstar_slice);
    alpha_eq_at_Lstar_slice = alpha_eq(:,:,nearest_Lstar_index);alpha_eq_at_Lstar_slice = alpha_eq_at_Lstar_slice(:);
    energy_at_Lstar_slice = energy(:,:,nearest_Lstar_index); energy_at_Lstar_slice = energy_at_Lstar_slice(:);
    flux_at_Lstar_slice = flux(:,:,nearest_Lstar_index); flux_at_Lstar_slice = flux_at_Lstar_slice(:);
    
    %remove flux<0 and to avoid log(0) errors, set flux = 0 to a very small value in comparison to typical flux values
    alpha_eq_at_Lstar_slice(flux_at_Lstar_slice<=0) = [];
    energy_at_Lstar_slice(flux_at_Lstar_slice<=0) = [];
    flux_at_Lstar_slice(flux_at_Lstar_slice<=0) = [];
    %flux_at_Lstar_slice(flux_at_Lstar_slice==0) = 10^(-40);
    logflux_at_Lstar_slice = log10(flux_at_Lstar_slice);
    
    %normalize before interpolation
    std_alpha = std(alpha_eq_at_Lstar_slice(:));
    std_energy = std(energy_at_Lstar_slice(:));
    mean_alpha = mean(alpha_eq_at_Lstar_slice(:));
    mean_energy = mean(energy_at_Lstar_slice(:));
    %normalize data
    alpha_eq_at_Lstar_slice = (alpha_eq_at_Lstar_slice-mean_alpha)/std_alpha;
    energy_at_Lstar_slice = (energy_at_Lstar_slice-mean_energy)/std_energy;
    %normalize interpolation points
    alpha_eq_interp_at_Lstar_slice = (alpha_eq_interp_at_Lstar_slice-mean_alpha)/std_alpha;
    energy_interp_at_Lstar_slice = (energy_interp_at_Lstar_slice-mean_energy)/std_energy;
    %added points used to avoid extrapolation
    alpha_extra_points = ([0;90]-mean_alpha)./std_alpha;
    energy_extra_points = ([4;4]-mean_energy)./std_energy;
    logflux_extra_points = [-41;-41];
    %interpolate (added points (flux = 0 at alpha = 0, 90 at E = 4MeV) to interpolate rather than extrapolate (the issue was that some extrapolation exploded and returned as Inf))
    warning('off')
    F = scatteredInterpolant([alpha_eq_at_Lstar_slice;alpha_extra_points],[energy_at_Lstar_slice;energy_extra_points],[logflux_at_Lstar_slice;logflux_extra_points],'linear','none');
    warning('on')
    logflux_interp_at_Lstar_slice = F(alpha_eq_interp_at_Lstar_slice,energy_interp_at_Lstar_slice);
    flux_interp_at_Lstar_slice = 10.^logflux_interp_at_Lstar_slice;
    %undo normalization
    alpha_eq_at_Lstar_slice = (alpha_eq_at_Lstar_slice*std_alpha)+mean_alpha;
    energy_at_Lstar_slice = (energy_at_Lstar_slice*std_energy)+mean_energy;
    alpha_eq_interp_at_Lstar_slice = (alpha_eq_interp_at_Lstar_slice*std_alpha)+mean_alpha;
    energy_interp_at_Lstar_slice = (energy_interp_at_Lstar_slice*std_energy)+mean_energy;
    
    if max(max(max(isinf(flux_interp_at_Lstar_slice)))) ~= 0
        disp('some flux components were interpolated as Inf!')
        disp(Lstar_slice)
    end
    
    %if not in the convex hull, interpolation returns NaN (change this to flux = 0)
    flux_interp_at_Lstar_slice(isnan(flux_interp_at_Lstar_slice)) = 0;
    
    %convert the interpolated flux to phase space density (PSD = flux/p^2)
    PSD_interp_at_Lstar_slice = flux_interp_at_Lstar_slice.*const.c^2./(energy_interp_at_Lstar_slice.^2+2*energy_interp_at_Lstar_slice*const.E_0);
    
    %interpolate weighting function and compute integral argument for integration
    W = interp1(rad2deg(alpha_eq_for_W_interp),W_for_W_interp,alpha_eq_interp_at_Lstar_slice);
    arg = abs(16.*(pi^2).*(const.R_E.^3).*(Lstar_slice.^2).*(const.E_0./(const.c.^3)).*((energy_interp_at_Lstar_slice./const.E_0) + 1).*sqrt(energy_interp_at_Lstar_slice.^2 + (2.*energy_interp_at_Lstar_slice.*const.E_0)).*W);
    
    %integrate over the pitch angle to obtain the differential content (TRBEC as a function of E and Lstar)
    diff_content(i,:) = trapz(alpha_eq_interp_at_Lstar_slice(:,1),arg.*PSD_interp_at_Lstar_slice,1);
    
    % %     %plots of the initial data and interpolated points, only shows first iteration to show that it is working as intended
    % %     if plot_toggle_interp == 1
    % %         figure()
    % %         scatter([alpha_eq_at_Lstar_slice;(alpha_extra_points*std_alpha+mean_alpha)],[energy_at_Lstar_slice;(energy_extra_points*std_energy+mean_energy)],20,[flux_at_Lstar_slice;10.^(logflux_extra_points)],'filled')
    % %         title(strcat('Data (L=',num2str(Lstar_slice),')'))
    % %         xlim([0 95])
    % %         ylim([0 1])
    % %         caxis([10^4 10^8])
    % %         xlabel('Pitch Angle (\circ)')
    % %         ylabel('Energy (MeV)')
    % %         set(gca,'colorscale','log')
    % %         cb = colorbar;
    % %         cb.Label.String = 'Flux';
    % %
    % %         figure()
    % %         scatter(alpha_eq_interp_at_Lstar_slice(:),energy_interp_at_Lstar_slice(:),20,flux_interp_at_Lstar_slice(:),'filled')
    % %         hold on
    % %         plot([alpha_eq_at_Lstar_slice;(alpha_extra_points*std_alpha+mean_alpha)],[energy_at_Lstar_slice;(energy_extra_points*std_energy+mean_energy)],'ko')
    % %         title(strcat('Interpolated (L=',num2str(Lstar_slice),')'))
    % %         xlim([0 95])
    % %         ylim([0 1])
    % %         caxis([10^4 10^8])
    % %         xlabel('Pitch Angle (\circ)')
    % %         ylabel('Energy (MeV)')
    % %         set(gca,'colorscale','log')
    % %         cb = colorbar;
    % %         cb.Label.String = 'Flux';
    % %     end
end


% % if plot_toggle_diff == 1
% %     %plot after integral over alpha as a function of E and L*
% %     figure()
% %     imagesc(Lstar_interp,energy_interp,diff_content');
% %     ax = gca;
% %     %date string on plot
% %     text(2.6,2,datestr(epoch_select,'mm/dd HH:MM'),'Color','w','FontSize',14)
% %     set(gca,'YDir','normal')
% %     set(gca,'colorscale','log')
% %     cb = colorbar;
% %     cb.Label.String = 'Differential Content (1/MeV)';
% %     set(gca, 'YScale', 'log')
% %     yticks([0.1 1])%use only 10^-1 and 10^0 as yticks
% %     yticklabels({'10^{-1}','10^0'})
% %     ax.XColor = 'w';%changes tickmarks to white
% %     ax.YColor = 'w';
% %     %change ticks labels back to black
% %     for i = 1:size(ax.XTickLabel)
% %         ax.XTickLabel{i} = ['\color{black}' ax.XTickLabel{i}];
% %     end
% %     for i = 1:size(ax.YTickLabel)
% %         ax.YTickLabel{i} = ['\color{black}' ax.YTickLabel{i}];
% %     end
% %     xlabel('L*','Color','k')
% %     ylabel('Energy (MeV)','Color','k')
% %     ylim([0 3])
% %     caxis([10^26 10^30])
% %     differential_content_filename = strcat('frame_',num2str(epoch_select),'.png');
% %     output_size = 1600*[1.5 1];
% %     resolution = 400;
% %     set(gcf,'Units','inches','Position',[0 0 output_size/resolution]);
% %
% %     exportgraphics(gcf,differential_content_filename,'Resolution',resolution)
% % end


%use if want to save plot as .png (specifically just the after integration
%plot so must have plot_toggle disabled and remove if statements around the plot above
%imwrite(getframe(gcf).cdata, strcat('C:\Users\Jared Pitzel\Desktop\GIF\',num2str(orbit_number),'.png'));%save plot
%H = getframe(gca);
%imwrite(H.cdata, strcat('C:\Users\Jared Pitzel\Desktop\GIF\',num2str(orbit_number),'.png'));%save axis

%integral before intgral over L*, for comparisons to mu,K space
%[Lstar_interp',trapz(energy_interp,diff_content,2)]

save(strcat(DC_directory,'DC_E_PA_',satellite,datestr(epoch_select,'_yyyymmdd_HHMMSS'),'.mat'),'diff_content','energy_interp','Lstar_interp','epoch_select')

%integrate differential content over energy and Lstar to obtain TRBEC
%TRBEC = trapz(energy_interp,trapz(Lstar_interp,diff_content,1),2);
end
