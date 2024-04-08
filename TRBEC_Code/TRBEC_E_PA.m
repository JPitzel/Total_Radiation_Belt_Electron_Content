
%What is does:
%Uses MagEIS electron flux data to compute the TRBEC in E,alpha,L* coordiantes

format longE
warning('off') 


tic
%names of data directories (flux data as a input, differential content as an output)
flux_directory  = 'flux_data/';
DC_directory = 'DC_E_PA/';
mkdir(DC_directory);%create direcotry to save differential to if it doesn't already exist



%start and stop days (inclues full half orbits that occured within the period)
epoch_start = datenum(2019,01,01);
epoch_stop = datenum(2019,12,31);

%integration bounds
energy_min = 0.1;
energy_max = 0.3;
Lstar_min = 3.5;
Lstar_max = 5;

%compute the TRBEC over desired time interval for given integration bounds
[TRBEC_rbspa,epoch_rbspa] = main(flux_directory,DC_directory,epoch_start,epoch_stop,'rbspa',energy_min,energy_max,Lstar_min,Lstar_max);
[TRBEC_rbspb,epoch_rbspb] = main(flux_directory,DC_directory,epoch_start,epoch_stop,'rbspb',energy_min,energy_max,Lstar_min,Lstar_max);
toc

%%
%get rid of TRBEC = 0
epoch_rbspa(TRBEC_rbspa == 0) = [];epoch_rbspb(TRBEC_rbspb == 0) = [];
TRBEC_rbspa(TRBEC_rbspa == 0) = [];TRBEC_rbspb(TRBEC_rbspb == 0) = [];

semilogy(epoch_rbspa,TRBEC_rbspa,'r')
hold on
semilogy(epoch_rbspb,TRBEC_rbspb,'b')


datetick('x','dd-mmm','keeplimits')
ylabel('TRBEC')
ax.XTickLabelMode = 'manual';
ax.XAxis.TickValues = datetime(2013,03,04)+calquarters(0:12);
ax.XAxis.TickLabelFormat = 'dd-mmm';
%epoch_start = datenum('2013-03-01 00:00','yyyy-mm-dd HH:MM');
%epoch_stop = datenum('2013-03-20 00:00','yyyy-mm-dd HH:MM');
xlim([epoch_start epoch_stop+1])
%ylim([10^23 10^29])
legend('RBSP A','RBSP B','Combined')


%add bounds to TRBEC plot
dim = [0.15, 0.6, 0.3, 0.3];
str = {[strcat('E =',32,num2str(energy_min),32,'to',32,num2str(energy_max),' MeV')],['\alpha_{eq} = \alpha_{eq,LC} to \pi/2'],[strcat('L* =',32,num2str(Lstar_min),32,'to',32,num2str(Lstar_max))]};
annotation('textbox',dim,'String',str,'FitBoxToText','on');



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    F = scatteredInterpolant([alpha_eq_at_Lstar_slice;alpha_extra_points],[energy_at_Lstar_slice;energy_extra_points],[logflux_at_Lstar_slice;logflux_extra_points],'linear','none');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%takes in data and selects out data for a given half orbit
function [flux_half_orbit,epoch_half_orbit,Lstar_half_orbit,B_sat_half_orbit,B_eq_half_orbit] = half_orbit_select(flux,epoch,Lstar,B_sat,B_eq,orbtimes,epoch_select)
%from the given time, select out the times at the beginning and end of the half orbit
[orbtimes_begin, orbtimes_begin_index]= max((orbtimes<=epoch_select).*orbtimes);
orbtimes_end = orbtimes(orbtimes_begin_index + 1);
%from the beginning and end times select out first and last index within the half orbit
%(NOTE: this assumes increasing values of time which is the case for epoch)
all_epoch_indices = find(epoch>=orbtimes_begin & epoch<orbtimes_end);%include start of orbit, exclude end (since it's the beginning of the next)
orbit_begin_index = all_epoch_indices(1);
orbit_end_index = all_epoch_indices(end);


%remove all data that isn't in the selected orbit based on the indices obtained above
flux_half_orbit = flux(:,:,orbit_begin_index:orbit_end_index);       
epoch_half_orbit = epoch(orbit_begin_index:orbit_end_index);     
Lstar_half_orbit = Lstar(orbit_begin_index:orbit_end_index);     
B_sat_half_orbit = B_sat(orbit_begin_index:orbit_end_index);
B_eq_half_orbit = B_eq(orbit_begin_index:orbit_end_index);

%remove all data that doesn't have an associated L* value
flux_half_orbit(:,:,Lstar_half_orbit<0) = [];
epoch_half_orbit(Lstar_half_orbit<0) = [];
B_sat_half_orbit(Lstar_half_orbit<0) = [];
B_eq_half_orbit(Lstar_half_orbit<0) = [];
Lstar_half_orbit(Lstar_half_orbit<0) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%computes the weight function as a function of alpha_eq, used for
%interpolation so no inputs required, output is in radians!
function [W,alpha_eq_mid] = weighting_function()

dalpha_eq = 1/100;%DO NOT TAKE TOO SMALL A STEP IN alpha_eq (causes issues when taking derivative of Y(y))
dtheta = 1/100000;%DO NOT TAKE A BIG STEP IN dtheta (causes issues when taking derivative of Y(y))

alpha_eq = (deg2rad(1):dalpha_eq:deg2rad(179))';

alpha_eq_length = length(alpha_eq);
y = sin(alpha_eq);

%need to go through each alpha_eq one at time, calculate theta_m
%and integrate over theta in order to obtain a y value,
%then once all y's have been obtained take derivative and are left with

%will need to compute T(y) for each alpha_eq bin, and for each time
%step since the loss cone changes over time
I_3 = zeros(alpha_eq_length,1);
for i = 1:alpha_eq_length
    theta_m = colatitude_m(alpha_eq(i));
    theta = theta_m:dtheta:pi/2;
    
    Y_y_theta = sin(theta).*sqrt(1+3.*(cos(theta)).^2).*((1-(y(i).^2.*sin(theta).^(-6).*sqrt(1+3*cos(theta).^2))).^(1/2));
    
    %take integral
    I_3(i) = (trapz(theta,Y_y_theta,2));
    
    %if Nan value is given (will occur at 90 deg) then replace with 0 (the actual value)
    if isnan(I_3(i))
       I_3(i) = 0; 
    end
end

%take derivative (taken from matlab documentation where derivative is approximated as diff()/dh)
I_3_derivative = diff(I_3)./dalpha_eq;
%derivative is estiamted at the mid points, not the original bin boundaries
alpha_eq_mid = (alpha_eq(2:end)+alpha_eq(1:end-1))/2;
%interpolate I_3 at the mid points (so it's the same size as the
%derivative)
I_3 = interp1(alpha_eq,I_3,alpha_eq_mid);

%calculate weighting function
W = abs(real(cos(alpha_eq_mid).*sin(alpha_eq_mid).*I_3 - sin(alpha_eq_mid).^2.*I_3_derivative));

%explicitly define W(0) = W(pi) = 0
alpha_eq_mid = [0;alpha_eq_mid;pi];
W = [0;W;0];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%given alpha_eq returns theta_m, input and output are in radians!
function theta_interp = colatitude_m(alpha_interp)
%since theta_m is symetric around 90, colat(alpha) = colat(180-alpha)
if alpha_interp > pi/2
    alpha_interp = pi - alpha_interp;
end
theta_m = 0:0.01:pi/2;
alpha_eq = asin(sqrt(((sin(theta_m).^6))./(sqrt(3*(cos(theta_m)).^2 + 1))));
theta_interp = interp1(alpha_eq,theta_m,alpha_interp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%computes half orbits given position data
%NOTE: marks times independently of units for position (some of the MagEIS data is given in R_E and others in km)
function orbtimes = get_orbtimes(filenames)
orbtimes = [];
for file_idx = 1:length(filenames)
    [data_orbtimes,~] = spdfcdfread(filenames{file_idx},'Variables',{'Epoch','Position'});
    epoch = data_orbtimes{1};
    position = data_orbtimes{2};
    distance = sqrt(sum(position.^2,2));
    
    %remove data that is too large
    epoch(distance>6*6371) = [];
    distance(distance>6*6371) = [];
    
    %if average is less than 10 R_E, convert to units of km
    %this fixes an issue where some versions of the MagEIS data (e.g. v8.4.0) use units of R_E instead of km as expected
    if mean(distance) < 10
        distance = distance*6371;
    end
    
    %smooth to avoid discontinuities that occur at apogee and perigee (i.e. there were instances where expected local min appeared as two local mins and one local max)
    distance = smooth(distance,21);
    
    %sort out local min and max that indicate satellite at perigee and apogee
    local_max_idx = islocalmax(distance);
    local_min_idx = islocalmin(distance);
    orbtimes_i = [epoch(local_max_idx);epoch(local_min_idx)];%combine and sort by occurance
    orbtimes_i = sortrows(orbtimes_i);

    orbtimes = cat(1,orbtimes,orbtimes_i);

    %plot distance vs time for troubleshooting
    %plot(epoch,distance)
    %hold on
end

%once the orbtimes are determined, remove orbits that are too short (some still occur even with smoothing done above)
orbperiod = diff(orbtimes);
idxs = find(orbperiod<0.95*mean(orbperiod));%find indices where difference between orbtimes is less than 95 percent of mean half orbit period
orbtimes(idxs+1) = [];%remove the identified orbtimes that were too short
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%given directory, satellite, and start and stop dates, outputs a list of
%files according to start/stop times and the epoch of each file as datenum
%satellite: either 'rbspa' or 'rbspb'
%assumes files are flux data, formatted as 'rbspx_rel04_ect-mageis-L3_yyyymmdd_vX.X.X.cdf'
function [selected_filenames,selected_file_epoch] = select_flux_filenames(flux_directory,satellite,start_epoch,end_epoch)

%for flux data, select out date from filename
%default formatting: 'rbspx_rel04_ect-mageis-L3_yyyymmdd_vX.X.X.cdf'
full_filenames = filelist(strcat(flux_directory,satellite,'*.*'));
%convert from string array (double quotes), to character vector (single quotes), necessary to strip date out of filenames
full_filenames_characters = char(full_filenames);
year = str2num(full_filenames_characters(:,27:30));
month = str2num(full_filenames_characters(:,31:32));
day = str2num(full_filenames_characters(:,33:34));
file_epoch = datenum(year,month,day);

%based on given dates, determines which files should be used and puts them in a list
selected_epoch_indices = find(file_epoch>=start_epoch & file_epoch<=end_epoch);%includes all full half orbits that occured within the time period
selected_filenames = full_filenames(selected_epoch_indices);
selected_file_epoch = file_epoch(selected_epoch_indices);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%given a directory (which can specify the type of file) returns a list of
%the files names in that folder (extracts the names from the "dir" function)
%example, filelist('data/PSD_rbspa*.*') will load all the files that start with 'PSD_rbspa' in the data folder
function output = filelist(directory)
files = dir(directory);
number_files = length(files);

file_names = strings(number_files,1);%create empty list
for i = 1:number_files
   file_names(i) = files(i).name;
end
output = file_names;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%from the differential content, computes the TRBEC
function [TRBEC, epoch] = TRBEC_energy_pitch_angle(DC_directory,satellite,selected_half_orbit,energy_min,energy_max,Lstar_min,Lstar_max)

load(strcat(DC_directory,'DC_E_PA_',satellite,'_',datestr(selected_half_orbit,'yyyymmdd_HHMMSS'),'.mat'));

if diff_content == 0
    TRBEC = 0;
    epoch = epoch_select;
else
%reduce data based on bounds
selected_energy_indices = find(energy_interp>=energy_min & energy_interp<=energy_max);
selected_Lstar_indices = find(Lstar_interp>=Lstar_min & Lstar_interp<=Lstar_max);

%compute the TRBEC using new bounds
TRBEC = trapz(energy_interp(selected_energy_indices),trapz(Lstar_interp(selected_Lstar_indices),diff_content(selected_Lstar_indices,selected_energy_indices),1),2);
epoch = epoch_select;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TRBEC,epoch] = main(flux_directory,DC_directory,start_epoch,stop_epoch,satellite,energy_min,energy_max,Lstar_min,Lstar_max)
%given start and stop times, determine filenames of flux data to import
[selected_flux_filenames,selected_flux_filename_epoch] = select_flux_filenames(flux_directory,satellite,start_epoch,stop_epoch);

%determine half orbit epochs
orbtimes = get_orbtimes(strcat(flux_directory,selected_flux_filenames));
TRBEC = []; epoch = [];

%for each half orbit compute the differential content and the TRBEC
for selected_half_orbit_idx = 1:length(orbtimes)-1%go through all complete half orbits
    selected_half_orbit = orbtimes(selected_half_orbit_idx);
    next_selected_half_orbit = orbtimes(selected_half_orbit_idx+1);
    
    %import only the necessary files to compute the DC/TRBEC for a single half orbit
    %(either 1 or 2 depending on if the half orbit spans multiple days)
    selected_file_idx = find(selected_flux_filename_epoch == floor(selected_half_orbit));%find filename index for half orbit start epoch
    %if next orbit is on the same day, only one file is need, otherwise the next day is also required to have data for the full half orbit
    if floor(selected_half_orbit) == floor(next_selected_half_orbit)
        filename_flux_i = selected_flux_filenames(selected_file_idx);
    else
        filename_flux_i = selected_flux_filenames(selected_file_idx:selected_file_idx+1);
    end
    
    %try loading the differential content to compute the TRBEC from it,
    %if a differential content file doesn't exist then it is generated
     try%try to compute TRBEC using DC file
         TRBEC = [TRBEC;TRBEC_energy_pitch_angle(DC_directory,satellite,selected_half_orbit,energy_min,energy_max,Lstar_min,Lstar_max)];
         epoch = [epoch;selected_half_orbit];
      catch
          try%if DC file doesn't already exist, generate it and compute the TRBEC from the DC file
             differential_content_energy_pitch_angle(selected_half_orbit,satellite,orbtimes,flux_directory,DC_directory,filename_flux_i);
             TRBEC = [TRBEC;TRBEC_energy_pitch_angle(DC_directory,satellite,selected_half_orbit,energy_min,energy_max,Lstar_min,Lstar_max)];
             epoch = [epoch;selected_half_orbit];
             disp(strcat('CREATED: ',DC_directory,'DC_E_A_',satellite,datestr(selected_half_orbit,'_yyyymmdd_HHMMSS'),'.mat'));
          catch%if error occurs for DC, then flag the half orbit
             disp('--------------------------------------------')
             disp(strcat('TRBEC_',satellite,', error at: ',datestr(selected_half_orbit,'yyyy-mm-dd HH:MM:SS')))
            
             %if error is thrown save the TRBEC as zero so it isn't computed every time
             diff_content = 0;
             epoch_select = selected_half_orbit;
             save(strcat(DC_directory,'DC_E_PA_',satellite,datestr(selected_half_orbit,'_yyyymmdd_HHMMSS'),'.mat'),'diff_content','epoch_select')
             disp(strcat('CREATED: ',DC_directory,'DC_E_A_',satellite,datestr(selected_half_orbit,'_yyyymmdd_HHMMSS'),'.mat'));
          end
     end
end
end
