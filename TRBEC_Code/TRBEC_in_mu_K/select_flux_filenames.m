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