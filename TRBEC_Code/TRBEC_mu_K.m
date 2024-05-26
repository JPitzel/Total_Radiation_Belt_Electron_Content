%Creation Date: February 25th, 2021
%What it does:
%
format longE
warning('off')


tic
%names of data directories (flux and PSD data as a inputs, differential content as an output)
PSD_directory  = 'PSD_data/';
flux_directory = 'flux_data/';
DC_directory = 'test_mu_K/';
mkdir(DC_directory);%create direcotry to save DC to if it doesn't already exist


%start and stop days (inclues full half orbits that occured within the period)
epoch_start = datenum(2013,03,01);
epoch_stop = datenum(2013,03,15);

%integration bounds
mu_min = 50;
mu_max = 200;
Lstar_min = 3.5;
Lstar_max = 5;

%compute the TRBEC over desired time interval forc given integration bounds
[TRBEC_rbspa,epoch_rbspa] = main(PSD_directory,flux_directory,DC_directory,epoch_start,epoch_stop,'rbspa',mu_min,mu_max,Lstar_min,Lstar_max);
[TRBEC_rbspb,epoch_rbspb] = main(PSD_directory,flux_directory,DC_directory,epoch_start,epoch_stop,'rbspb',mu_min,mu_max,Lstar_min,Lstar_max);
toc

%%
%remove TRBEC = 0
epoch_rbspa(TRBEC_rbspa == 0) = [];epoch_rbspb(TRBEC_rbspb == 0) = [];
TRBEC_rbspa(TRBEC_rbspa == 0) = [];TRBEC_rbspb(TRBEC_rbspb == 0) = [];
%remove outliers ('mean' defines outliers as points that are three standard deviations from the mean).
[TRBEC_rbspa_cleaned,rbspa_cleaned_indices] = rmoutliers(TRBEC_rbspa);[TRBEC_rbspb_cleaned,rbspb_cleaned_indices] = rmoutliers(TRBEC_rbspb);
epoch_rbspa_cleaned = epoch_rbspa(rbspa_cleaned_indices == 0); epoch_rbspb_cleaned = epoch_rbspb(rbspb_cleaned_indices == 0); 


semilogy(epoch_rbspa_cleaned,TRBEC_rbspa_cleaned,'r-')
hold on
semilogy(epoch_rbspb_cleaned,TRBEC_rbspb_cleaned,'b-')


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
str = {[strcat('\mu =',32,num2str(mu_min),32,'to',32,num2str(mu_max),' MeV')],['\alpha_{eq} = \alpha_{eq,LC} to \pi/2'],[strcat('L* =',32,num2str(Lstar_min),32,'to',32,num2str(Lstar_max))]};
annotation('textbox',dim,'String',str,'FitBoxToText','on');


% % output_size = 1500*[3.5 1];
% % resolution = 500;
% % set(gcf,'Units','inches','Position',[0 0 output_size/resolution]);
% % 
% % filename_export = strcat('adiabatic_TRBEC_',num2str(mu_min),'_to_',num2str(mu_max),'.png');
% % exportgraphics(gcf,filename_export,'Resolution',resolution)
% % xticks(datenum(datetime('21-Feb-2013'):caldays(2):datetime('21-Mar-2013')))
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%given half orbit and filenames to import data from, outputs the differential content
%as a function of first and third adiabatic invariants (mu,L*)
function differential_content_mu_K(epoch_select,satellite,orbtimes,PSD_directory,flux_directory,DC_directory,filename_PSD,filename_flux)

%constants
const.R_E = 6371; %km, Earth radius
const.E_0 = 0.511; %MeV, electron rest energy
const.B_0 = 0.311653; %G
const.mu_E = const.B_0*(const.R_E*10^5)^3; %cm^3 G, magnetic moment of Earth

%import PSD data from given file list
PSD = []; epoch = []; mu = []; K = []; Lstar = []; energy = []; alpha = [];
for PSD_file_idx = 1:length(filename_PSD)
    [PSD_i,epoch_i,mu_i,K_i,Lstar_i,energy_i,alpha_i,datainfo] = import_PSD(strcat(PSD_directory,filename_PSD(PSD_file_idx)));
    PSD = cat(3,PSD,PSD_i);
    epoch = cat(1,epoch,epoch_i);
    mu = cat(3,mu,mu_i);
    K = cat(3,K,K_i);
    Lstar = cat(3,Lstar,Lstar_i);
    energy = cat(3,energy,energy_i);
    alpha = cat(3,alpha,alpha_i);
end

%selects out the data for the half orbit that occured during the selected time
[PSD,epoch,mu,K,Lstar,energy,alpha] = half_orbit_select(PSD,epoch,mu,K,Lstar,energy,alpha,orbtimes,epoch_select);
[K_length,mu_length,Lstar_length] = size(PSD);

%import B field data from flux data (its interpolated to be synced with PSD data)
%(note: if moved from this spot in the code must change alpha = 90 index in the get_B_fields function)
[Bsat,Beq] = get_B_fields(strcat(flux_directory,filename_flux),epoch,Lstar);

%add nearest point estimate for the loss cone at small mu
[PSD,mu,K,Lstar,alpha,energy] = K_loss_cone_small_mu(PSD,mu,K,Lstar,alpha,energy,Bsat,const,1);%using K(alpha = 10)
[PSD,mu,K,Lstar,alpha,energy] = K_loss_cone_small_mu(PSD,mu,K,Lstar,alpha,energy,Bsat,const,length(alpha(:,1,1))-1);%using K(alpha = 170)

%add PSD = 0 at the loss cone at large mu
[PSD,mu,K,Lstar,alpha,energy] = K_loss_cone_large_mu(PSD,mu,K,Lstar,alpha,energy);

%add in K = 0 particles that cannot be seen by the detector
[PSD,mu,K,Lstar,alpha,energy] = equatorial_particles(PSD,mu,K,Lstar,alpha,energy,Bsat,Beq);

%sets K as a placeholder value when the PSD is a placeholder (these points are then ignored in the interpolation)
K(PSD == 0) = -10^31;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%select bounds and resolution for interpolation (NOTE: K is from 0 to K_LC)
mu_min = 50;
mu_max = 2000;
mu_spacing = 10;
Lstar_min = 3.5;
Lstar_max = 5;
Lstar_spacing = 0.1;
K_interp_length = 200;
mu_interp_length = round((mu_max-mu_min)/mu_spacing + 1);
Lstar_interp_length = length((Lstar_min:Lstar_spacing:Lstar_max));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%interpolates over the region of interest to compute the differential content
%and saves it to a seperate file
%use log10(PSD) for interpolation, but first change where PSD = 0 to avoid log(0) issues
PSD(PSD == 0) = 10^-40;
logPSD = log10(PSD);

%finds mean Lstar at each time, used to find closest value later on
[~,~,Lstar_length] = size(Lstar);
Lstar_mean = zeros(Lstar_length,1);
for i = 1:Lstar_length
    Lstar_mean(i) = mean(mean(Lstar(:,:,i)));
end

Lstar_limits = linspace(Lstar_min,Lstar_max,Lstar_interp_length);
mu_limits = logspace(log10(mu_min),log10(mu_max),mu_interp_length);


diff_content = ones(mu_interp_length,Lstar_interp_length);

for i = 1:length(Lstar_limits)
    
    Lstar_slice = Lstar_limits(i);
    [~,nearest_Lstar_index] = min(abs(mean(mean(Lstar-Lstar_slice))));
    
    K_at_Lstar_slice = K(:,:,nearest_Lstar_index); K_at_Lstar_slice = K_at_Lstar_slice(:);
    mu_at_Lstar_slice = mu(:,:,nearest_Lstar_index); mu_at_Lstar_slice = mu_at_Lstar_slice(:);
    logPSD_at_Lstar_slice = logPSD(:,:,nearest_Lstar_index); logPSD_at_Lstar_slice = logPSD_at_Lstar_slice(:);
    
    mu_at_Lstar_slice(K_at_Lstar_slice<0) = [];
    logPSD_at_Lstar_slice(K_at_Lstar_slice<0) = [];
    K_at_Lstar_slice(K_at_Lstar_slice<0) = [];
    K_at_Lstar_slice(mu_at_Lstar_slice<0) = [];
    logPSD_at_Lstar_slice( mu_at_Lstar_slice<0) = [];
    mu_at_Lstar_slice(mu_at_Lstar_slice<0) = [];
    
% %     plot_toggle = 0;
% %     if plot_toggle == 1
% %         figure()
% %         scatter(K_at_Lstar_slice,mu_at_Lstar_slice,20,10.^logPSD_at_Lstar_slice,'filled')
% %         title(strcat('interpolated','(L* \approx',num2str(Lstar_slice),')'))
% %         set(gca,'colorscale','log')
% %         %xlim([0 2.11020679*Lstar_slice-2.26735643])
% %         xlim([0 0.5])
% %         ylim([50 200])
% %         xlabel('K (R_E \surdG)','Interpreter','tex')
% %         ylabel('\mu (MeV/G)')
% %         cb = colorbar;
% %         cb.Label.String = 'PSD [(c/(cm MeV))^3]';
% %         caxis([10^-4 10^-2])
% %     end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %sets values of K interpolation
    K_upper_limit = 2.11020679*Lstar_limits(i)-2.26735643; %K_LC(Lstar) function
    K_limits = logspace(-41,log10(K_upper_limit),K_interp_length);%K range is from 0 (or close it it on a log scale) to K_LC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %calculate standard deviations and means of the data for use in normalization and reverse process later on
    std_K = std(K_at_Lstar_slice(:));
    std_mu = std(mu_at_Lstar_slice(:));
    mean_K = mean(K_at_Lstar_slice(:));
    mean_mu = mean(mu_at_Lstar_slice(:));
    %normalize K and mu data
    K_at_Lstar_slice = (K_at_Lstar_slice - mean_K)/std_K;
    mu_at_Lstar_slice = (mu_at_Lstar_slice - mean_mu)/std_mu;
    %normalize K and mu interpolation points
    [K_interp,mu_interp] = meshgrid(K_limits,mu_limits);
    K_interp = (K_interp - mean_K)/std_K;
    mu_interp = (mu_interp - mean_mu)/std_mu;
    
    %interpolate
    F = scatteredInterpolant(K_at_Lstar_slice(:),mu_at_Lstar_slice(:),logPSD_at_Lstar_slice(:),'linear','linear');
    logPSD_interp = F(K_interp,mu_interp);
    
    %flag if interpolation isn't as expected
    if sum(sum(sum(isnan(logPSD_interp)))) ~= 0
        disp('some PSD components were interpolated as NaN!')
        disp(Lstar_slice)
    end
    
    %unnormalize interpolated points and recover PSD from log10(PSD)
    K_interp = (K_interp*std_K) + mean_K;
    mu_interp = (mu_interp*std_mu) + mean_mu;
    PSD_interp = 10.^logPSD_interp;
    
    
    %another flag if the interpolation isn't as expected
    if max(max(max(isinf(PSD_interp)))) ~= 0
        disp('some PSD components were interpolated as Inf!')
        disp(Lstar_slice)
    end
    
    

    %calculates the TRBEC integral argument (neglecting constants)
    TRBEC_int_arg = sqrt(mu_interp).*PSD_interp./(Lstar_slice.^2);
    
    %compute differential content by integrating over K
    diff_content(:,i) = 8.*sqrt(2).*(pi^2).*(const.E_0^(3/2)).*const.mu_E.*trapz(K_limits,TRBEC_int_arg,2);
    

    % %     %contour plots of integral argument vs mu and K at varying Lstar
    % %         if plot_toggle == 1
    % %             figure()
    % %             scatter(K_interp(:),mu_interp(:),20,PSD_interp(:),'filled')
    % %             title(strcat('interpolated','(L* \approx',num2str(Lstar_slice),')'))
    % %             set(gca,'colorscale','log')
    % %             %xlim([0 2.11020679*Lstar_slice-2.26735643])
    % %             xlim([0 0.5])
    % %             ylim([50 200])
    % %             xlabel('K (R_E \surdG)','Interpreter','tex')
    % %             ylabel('\mu (MeV/G)')
    % %             cb = colorbar;
    % %             cb.Label.String = 'PSD [(c/(cm MeV))^3]';
    % %             caxis([10^-4 10^-2])
    % %             hold on
    % %             contour(K_interp,mu_interp,PSD_interp,(10^-4:8*10^-4:10^-2),'k')
    % %         end
end

%saves the differential content to a seperate file
mu_limits = mu_limits';Lstar_limits = Lstar_limits';%invert these vectors now, unsure if changing it earlier would affect code
save(strcat(DC_directory,'DC_mu_K_',satellite,datestr(epoch_select,'_yyyymmdd_HHMMSS'),'.mat'),'diff_content','mu_limits','Lstar_limits','epoch_select')

%plots after integral over K, integrated integral arguments vs mu and Lstar

% % if plot_toggle == 1
% %     figure()
% %     imagesc(Lstar_limits,mu_limits,diff_content);
% %     title('after K integral')
% %     set(gca,'YDir','normal')
% %     set(gca,'colorscale','log')
% %     cb = colorbar;
% %     cb.Label.String = '';
% %     xlabel('L*')
% %     ylabel('\mu (MeV/G)')
% %     cb.Label.String = 'PSD Integrated Over K';
% % end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%imports data, splits into variables and changes index order to match 
%(dim 1 = # pitch angle bins, dim 2 = # energy bins, dim 3 = # epoch bins)
%and adds dimensions where applicable (for instance K is presented in 2D grid, but it's converted to a 3D one)
function [PSD,epoch,mu,K,Lstar,energy,alpha,datainfo] = import_PSD(filename)
[data,datainfo] = spdfcdfread(filename,'Variable',{'PSD','Epoch','mu','I','Lstar','EL','PA'});

%NOTE:on rbspgway I is K
%energy: even though it says keV in the metadata, it's MeV (lowest bin on MagEIS is known to be 38keV and yet is listed as 3.8*10^-2)
PSD = data{1};epoch = data{2};mu = data{3};K = data{4};Lstar = data{5};energy = data{6};alpha = data{7};
[number_alpha_bins,number_energy_bins,~] = size(PSD);

%for ease of use put all variables in a meshgrid except epoch (NOTE: PSD, mu are in a 3D grid by default)
PSD = double(PSD);
K = repmat(K, [1,1,number_energy_bins]); K = permute(K,[2,3,1]);
mu = double(mu);
Lstar = repmat(Lstar,[1,1,number_energy_bins]); Lstar = permute(Lstar,[2,3,1]);
energy = repmat(energy,[1,1,number_alpha_bins]); energy = permute(energy,[3,2,1]);
alpha = repmat(alpha,[1,1,number_energy_bins]); alpha = permute(alpha,[2,3,1]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%given directory, satellite, and start and stop dates, outputs a list of
%files according to start/stop times and the epoch of each file as datenum
%satellite: either 'rbspa' or 'rbspb'
%assumes files are flux data, formatted as 'PSD_rbspx_mageis_yyyyddd.cdf'
function [selected_filenames,selected_file_epoch] = select_PSD_filenames(data_directory,satellite,start_time,end_time)
%get dates from PSD filenames
%for PSD data, select out date from filename
%default formatting: 'PSD_rbspx_mageis_yyyyddd.cdf'
full_filenames = filelist(strcat(data_directory,'PSD_',satellite,'*.*'));
%convert from string array (double quotes), to character vector (single quotes), necessary to strip date out of filenames
full_filenames_characters = char(full_filenames);
year = str2num(full_filenames_characters(:,18:21));
doy = str2num(full_filenames_characters(:,22:24));
file_epoch = datenum(year, ones(size(year)), doy);

%based on given dates, determines which files should be used and puts them in a list
selected_epoch_indices = find(file_epoch>=start_time & file_epoch<=end_time);%includes all full half orbits that occured within the time period
selected_filenames = full_filenames(selected_epoch_indices);
selected_file_epoch = file_epoch(selected_epoch_indices);
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

%Takes in a given time and orbit data and return data for the half orbit that occured during the given time
function [PSD_half_orbit,epoch_half_orbit,mu_half_orbit,K_half_orbit,Lstar_half_orbit,energy_half_orbit,alpha_half_orbit] = half_orbit_select(PSD,epoch,mu,K,Lstar,energy,alpha,orbtimes,epoch_select)
%from the given time, select out the times that denote the beginning and end of the half orbit
[orbtimes_begin, orbtimes_begin_index]= max((orbtimes<=epoch_select).*orbtimes);
orbtimes_end = orbtimes(orbtimes_begin_index + 1);
%from the beginning and end times select out first and last index within the half orbit
%(NOTE: this assumes increasing values of time which is the case for epoch)
all_epoch_indices = find(epoch>=orbtimes_begin & epoch<orbtimes_end);%include start of orbit, exclude end
orbit_begin_index = all_epoch_indices(1);
orbit_end_index = all_epoch_indices(end);

%remove all data that isn't in the selected orbit based on the indices obtained above
PSD_half_orbit = PSD(:,:,orbit_begin_index:orbit_end_index);       
epoch_half_orbit = epoch(orbit_begin_index:orbit_end_index);     
mu_half_orbit = mu(:,:,orbit_begin_index:orbit_end_index);         
K_half_orbit = K(:,:,orbit_begin_index:orbit_end_index);
Lstar_half_orbit = Lstar(:,:,orbit_begin_index:orbit_end_index);     
energy_half_orbit = energy(:,:,orbit_begin_index:orbit_end_index);
alpha_half_orbit = alpha(:,:,orbit_begin_index:orbit_end_index);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%adds in points for the loss cone (alpha_LC converted to K_LC) according to steps outlined in overleaf documentation
function [PSD,mu,K,Lstar,alpha,energy] = K_loss_cone_small_mu(PSD,mu,K,Lstar,alpha,energy,Bsat,const,alpha_idx)
[~,mu_length,Lstar_length] = size(PSD);
%step 1
h = 100; %km, choosen altitude at which particles precipitate
alpha_LC = asind((1 + ((1.5).*(h./const.R_E)))./(Lstar(alpha_idx,:,:).^(3/2).*(4-((3./Lstar(alpha_idx,:,:)).*(1+(h./const.R_E)))).^(0.25))); %equation from JF's work
K_LC = 2.11020679*Lstar(alpha_idx,:,:)-2.26735643; %equation from Chris
energy_LC = energy(alpha_idx,:,:);

%step 2
mu_LC = zeros(1,mu_length,Lstar_length);
for i = 1:Lstar_length
    mu_LC(1,:,i) = ((energy_LC(1,:,i).^2+2*const.E_0.*energy_LC(1,:,i)).*(sind(alpha_LC(1,:,i)).^2))./(2*const.E_0.*Bsat(i));
end

Lstar_LC = zeros(1,mu_length,Lstar_length);
for mu_idx = 1:mu_length
    for Lstar_idx = 1:Lstar_length
        if isreal(alpha_LC(1,mu_idx,Lstar_idx))
            Lstar_LC(1,mu_idx,Lstar_idx) = interp1(alpha(:,mu_idx,Lstar_idx),Lstar(:,mu_idx,Lstar_idx),alpha_LC(1,mu_idx,Lstar_idx),'linear','extrap');
        else
            alpha_LC(1,mu_idx,Lstar_idx) = -10^31;
            Lstar_LC(1,mu_idx,Lstar_idx) = -10^31;
        end
    end
end
PSD_LC = zeros(1,mu_length,Lstar_length);
for mu_idx = 1:mu_length
    for Lstar_idx = 1:Lstar_length
        if alpha_LC(1,mu_idx,Lstar_idx) > 0
            PSD_LC(1,mu_idx,Lstar_idx) = PSD(alpha_idx,mu_idx,Lstar_idx);
        else
            PSD_LC(1,mu_idx,Lstar_idx) = 0;
        end
    end
end
%concatenate the loss cone with the RBSP data
PSD = cat(1,PSD,PSD_LC);
mu = cat(1,mu,mu_LC);
K = cat(1,K,K_LC);
Lstar = cat(1,Lstar,Lstar_LC);
alpha = cat(1,alpha,alpha_LC);
energy = cat(1,energy,energy_LC);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%adds in PSD = 0 at the loss cone for large mu
function [PSD,mu,K,Lstar,alpha,energy] = K_loss_cone_large_mu(PSD,mu,K,Lstar,alpha,energy)
[~,mu_length,Lstar_length] = size(PSD);
%add in K_LC for high values of mu (make interpolatation easier and hence the integral from 0 to K_LC
PSD = cat(1,PSD,10^-40*ones(1,mu_length,Lstar_length));%use small number for PSD and not zero because PSD = 0 is fill value of data
mu = cat(1,mu,3000*ones(1,mu_length,Lstar_length));
K = cat(1,K,K(end,:,:));
Lstar = cat(1,Lstar,Lstar(end,:,:));
alpha = cat(1,alpha,zeros(1,mu_length,Lstar_length));%values are meaningless, but need to something so dimensions are consistent between every variable
energy = cat(1,energy,zeros(1,mu_length,Lstar_length));%these are also meaningless
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%add in K = 0 particles that cannot be seen by the detector
%(NOTE: small K corresponds to pitch angles close to 90 so index 17 is alpha = 90 which is the smallest K)
function [PSD,mu,K,Lstar,alpha,energy] = equatorial_particles(PSD,mu,K,Lstar,alpha,energy,Bsat,Beq)
[~,mu_length,Lstar_length] = size(PSD);
K = cat(1,zeros(1,mu_length,Lstar_length),K);
mu = cat(1,mu(17,:,:),mu);
Lstar = cat(1,Lstar(17,:,:),Lstar);
PSD = cat(1,PSD(17,:,:),PSD);
%PSD = cat(1,zeros(1,mu_length,Lstar_length),PSD);%sets PSD = 0 for testing, switch with line above
alpha = cat(1,90*ones(1,mu_length,Lstar_length),alpha);%K = 0 corresponds to alpha = 90;
energy = cat(1,energy(17,:,:),energy);%not sure, but energies should be close to the energies at 90

%corrects all mu values at K = 0 from nearest point estimate to better
%approximation using curves of constant energy
for i = 1:Lstar_length
    mu(1,:,i) = Bsat(i)./Beq(i).*mu(1,:,i);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%given epochs and Lstar of given half orbit, imports and interpolates B field data
function [Bsat_interp,Beq_interp] = get_B_fields(flux_filenames,epoch,Lstar)

%import data
data = [];
for i = 1:length(flux_filenames)
    [data_i,~] = spdfcdfread(flux_filenames(i),'Variables',{'Epoch','L_star','B_Calc','B_Eq'});
    data_i = [data_i{1},data_i{2},data_i{3}*10^-5,data_i{4}*10^-5];%*10^-5 is to convert from nT to G
    data = [data;data_i];
end

%select out start and stop times of half orbit
start_half_orbit = epoch(1);
end_half_orbit = epoch(end);

%select out 90 deg pitch angle Lstars from PSD Lstar values (since Lstar from CDAweb is calculated using 90 deg particles)
%NOTE: based on the placement of the function is the main code index 17 is alpha = 90
Lstar_interp = squeeze(Lstar(17,1,:));%squeeze removes dimensions that are 1


%remove data where L_star, or B is less than zero and select times only within the half orbit
data(data(:,1)< start_half_orbit,:) = [];
data(data(:,1)>= end_half_orbit,:) = [];
data(data(:,2)<0,:) = [];
data(data(:,3)<0,:) = [];
data(data(:,4)<0,:) = [];

%seperate individual elements
Lstar_from_B_fields = data(:,2);
Bsat = data(:,3); %G
Beq = data(:,4); %G

%interpolate Bsat and Beq at the same points as the PSD data
Bsat_interp = interp1(Lstar_from_B_fields,Bsat,Lstar_interp);
Beq_interp = interp1(Lstar_from_B_fields,Beq,Lstar_interp);

Bsat_interp(isnan(Bsat_interp)) = -10^31;
Beq_interp(isnan(Beq_interp)) = -10^31;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%given filesnames of the differential content and integration bounds, computes the TRBEC
function [TRBEC,epoch] = TRBEC_adiabatic(DC_directory,satellite,selected_half_orbit,mu_min,mu_max,Lstar_min,Lstar_max)

load(strcat(DC_directory,'DC_mu_K_',satellite,'_',datestr(selected_half_orbit,'yyyymmdd_HHMMSS'),'.mat'))

%if diff_content is intentionally set to zero (scalar), then set TRBEC to zero
if diff_content == 0
TRBEC = 0;
epoch = epoch_select; 
else
%reduce data based on bounds
selected_mu_indices = find(mu_limits>=mu_min & mu_limits<=mu_max);
selected_Lstar_indices = find(Lstar_limits>=Lstar_min & Lstar_limits<=Lstar_max);

%compute the TRBEC using new bounds
TRBEC = trapz(mu_limits(selected_mu_indices),trapz(Lstar_limits(selected_Lstar_indices),diff_content(selected_mu_indices,selected_Lstar_indices),2),1);
epoch = epoch_select;
end
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

function [TRBEC, epoch] = main(PSD_directory,flux_directory,DC_directory,start_epoch,stop_epoch,satellite,mu_min,mu_max,Lstar_min,Lstar_max)
%given start and stop times, determine filenames of flux and PSD data that need to be imported
[selected_PSD_filenames,selected_PSD_filename_epoch] = select_PSD_filenames(PSD_directory,satellite,start_epoch,stop_epoch);
[selected_flux_filenames,~] = select_flux_filenames(flux_directory,satellite,start_epoch,stop_epoch);

%determine half orbit epochs
orbtimes = get_orbtimes(strcat(flux_directory,selected_flux_filenames));

TRBEC = []; epoch = [];
%for each half orbit compute the differential content and the TRBEC
for selected_half_orbit_idx = 1:length(orbtimes)-1%go through all complete half orbits
    selected_half_orbit = orbtimes(selected_half_orbit_idx);
    next_selected_half_orbit = orbtimes(selected_half_orbit_idx+1);
    
    %import only the necessary files to compute the DC/TRBEC for a single half orbit
    %(either 1 or 2 depending on if the half orbit spans multiple days)
    selected_file_idx = find(selected_PSD_filename_epoch == floor(selected_half_orbit));%find filename index for half orbit start epoch
    %if next orbit is on the same day, only one file is need, otherwise the next day is also required to have data for the full half orbit
    if floor(selected_half_orbit) == floor(next_selected_half_orbit)
        filename_PSD_i = selected_PSD_filenames(selected_file_idx);
        filename_flux_i = selected_flux_filenames(selected_file_idx);
    else
        filename_PSD_i = selected_PSD_filenames(selected_file_idx:selected_file_idx+1);
        filename_flux_i = selected_flux_filenames(selected_file_idx:selected_file_idx+1);
    end
    
    %try loading the differential content to compute the TRBEC from it,
    %if a differential content file doesn't exist then it is generated
    try%try to compute TRBEC using DC file
        TRBEC = [TRBEC;TRBEC_adiabatic(DC_directory,satellite,selected_half_orbit,mu_min,mu_max,Lstar_min,Lstar_max)];
        epoch = [epoch;selected_half_orbit];
    catch
        try%if DC file doesn't already exist, generate it and compute the TRBEC from the DC file
            differential_content_mu_K(selected_half_orbit,satellite,orbtimes,PSD_directory,flux_directory,DC_directory,filename_PSD_i,filename_flux_i)
            
            TRBEC = [TRBEC;TRBEC_adiabatic(DC_directory,satellite,selected_half_orbit,mu_min,mu_max,Lstar_min,Lstar_max)];
            epoch = [epoch;selected_half_orbit];
            disp(strcat('CREATED: ',DC_directory,'DC_mu_K_',satellite,datestr(selected_half_orbit,'_yyyymmdd_HHMMSS'),'.mat'));
        catch%if error occurs for DC, then flag the half orbit
            disp('---------------------------------------------------')
            disp(strcat('TRBEC_',satellite,', error at: ',datestr(selected_half_orbit,'yyyy-mm-dd HH:MM:SS')))
            
            %save file for TRBEC will be zero so don't compute orbits that give errors every time the code runs
            diff_content = 0;
            epoch_select = selected_half_orbit;
            save(strcat(DC_directory,'DC_mu_K_',satellite,datestr(epoch_select,'_yyyymmdd_HHMMSS'),'.mat'),'diff_content','epoch_select')
            disp(strcat('CREATED: ',DC_directory,'DC_mu_K_',satellite,datestr(selected_half_orbit,'_yyyymmdd_HHMMSS'),'.mat'));
  
        end
    end
    
              
end
end
