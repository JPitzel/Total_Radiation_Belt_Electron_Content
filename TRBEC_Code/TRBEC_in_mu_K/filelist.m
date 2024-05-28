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