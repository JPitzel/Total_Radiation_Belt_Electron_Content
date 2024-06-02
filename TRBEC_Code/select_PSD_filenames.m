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