
function [TRBEC, epoch] = TRBEC_mu_K_wrapper(PSD_directory,flux_directory,DC_directory,start_epoch,stop_epoch,satellite,mu_min,mu_max,Lstar_min,Lstar_max)
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
