
function [TRBEC,epoch] = TRBEC_energy_pitch_angle_wrapper(flux_directory,DC_directory,start_epoch,stop_epoch,satellite,energy_min,energy_max,Lstar_min,Lstar_max)
% Compute the TRBEC using (Energy, L*) bounds.
% Use the differential content file if it exists.
% Otherwise, generate the differential content file first.

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