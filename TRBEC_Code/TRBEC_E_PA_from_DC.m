%Creation Date: February 28th, 2023
%What is does:
%uses DC (output of TRBEC_E_PA.m) and plots the TRBEC given many sets of integration bounds

clear
clc

tic
directory = 'DC_E_PA\';
export_filename = 'E_PA_35_to_4.jpg';
legend_toggle = 1;
bounds_toggle = 1;
output_toggle = 1;

%start and stop days (inclues half orbits that started on the given dates)
epoch_start = datenum(2013,02,27);
epoch_stop = datenum(2013,03,20);

%for each set of integration bounds, combined data or keep seperate
separate_or_combine = 'sep';%'sep' to keep satellite data separate, 'comb' to combine data together (is also smoothed)
%for each set of integration bounds, single or multiple plots
individual_or_multiple = 'sgl'; %'sgl' for individual or 'mult' for multiple


colors = [0.000, 0.447, 0.741;
           0.850, 0.325, 0.098;
           0.929, 0.694, 0.125;	
           0.494, 0.184, 0.556;
           0.466, 0.674, 0.188;
           0.301, 0.745, 0.933;	
           0.635, 0.078, 0.184;];

%energy bounds [E_min,E_max]
E_bounds = [0.05,0.1;
            0.1, 0.3;
            0.3, 1.0];
%Lstar bounds [Lstar_min,Lstar_max]
Lstar_min = 3.5;
Lstar_max = 5;
Lstar_bounds = [Lstar_min,Lstar_max;
                Lstar_min,Lstar_max;
                Lstar_min,Lstar_max;];

main(directory,epoch_start,epoch_stop,separate_or_combine,individual_or_multiple,colors,E_bounds,Lstar_bounds)


fontsize = 20;
font_type = 'times';
set(gca,'fontname',font_type)

%axis labels
ylabel('TRBEC')
xlim([epoch_start epoch_stop]);
ax = gca; 
ax.FontSize = fontsize; 
set(gca,'yscale','log')%explicitly ensure plots is y-scale, issue where it would be linear even using semilogy!!!
xticks(datenum(datetime(datevec(epoch_start)):caldays(2):datetime(datevec(epoch_stop))))
datetick('x','dd-mmm','keeplimits','keepticks')
%ylim([10^25 10^29])

if legend_toggle == 1
    %legends (bounds for E)
    if strcmp(separate_or_combine,'sep') == 1
        legend_plots{1} = plot(nan,'k-');legend_labels{1} = "MagEIS A";
        legend_plots{2} = plot(nan,'k-.');legend_labels{2} = "MagEIS B";
        legend_length = 2;
    else
        legend_length = 0;
    end
    for i = (1:length(E_bounds))
        legend_plots{i+legend_length} = plot(nan,'Color',colors(i,:)) ;
        legend_labels{i+legend_length} = [num2str(E_bounds(i,1)) ' to ' num2str(E_bounds(i,2)) ' MeV'];
    end
    dim_legend = [0.35, 0.68, 0.2, 0.1];
    [~, hobj, ~, ~] = legend([legend_plots{:}],legend_labels{:},'location',dim_legend,'FontSize',fontsize);
    hl = findobj(hobj,'type','line');
    set(hl,'LineWidth',2);
end

if bounds_toggle == 1
%bounds for PA and L*
dim_bounds = [0.55, 0.6045, 0.3, 0.3];
str = {['\alpha_{eq} = \alpha_{eq,LC} to \pi/2'],[strcat('L* =',32,num2str(Lstar_bounds(1,1)),32,'to',32,num2str(Lstar_bounds(1,2)))]};
annotation('textbox',dim_bounds,'String',str,'FitBoxToText','on','FontSize',fontsize,'BackgroundColor','w','FaceAlpha',1,'FontName',font_type);
end

if output_toggle == 1
    %output
    resolution = 500;
    exportgraphics(gcf,export_filename,'Resolution',resolution)
end


toc


function main(directory,epoch_start,epoch_stop,separate_or_combine,individual_or_multiple,colors,E_bounds,Lstar_bounds)

[filenames_rbspa,~] = select_data_filenames(directory,'rbspa',epoch_start,epoch_stop);
[filenames_rbspb,~] = select_data_filenames(directory,'rbspb',epoch_start,epoch_stop);

[number_of_bounds,~] = size(E_bounds);

for j = 1:number_of_bounds

    [TRBEC_rbspa,epoch_rbspa] = TRBEC_E_PA(directory,filenames_rbspa,E_bounds(j,1),E_bounds(j,2),Lstar_bounds(j,1),Lstar_bounds(j,2));
    [TRBEC_rbspb,epoch_rbspb] = TRBEC_E_PA(directory,filenames_rbspb,E_bounds(j,1),E_bounds(j,2),Lstar_bounds(j,1),Lstar_bounds(j,2));
    
    %get rid of TRBEC = 0
    epoch_rbspa(TRBEC_rbspa == 0) = [];epoch_rbspb(TRBEC_rbspb == 0) = [];
    TRBEC_rbspa(TRBEC_rbspa == 0) = [];TRBEC_rbspb(TRBEC_rbspb == 0) = [];

    if strcmp(individual_or_multiple,'mult') == 1
       figure() 
    end
    
    if strcmp(separate_or_combine,'sep')
        
        semilogy(epoch_rbspa,TRBEC_rbspa,'Color',colors(j,:),'LineStyle','-','LineWidth',2)
        hold on
        semilogy(epoch_rbspb,TRBEC_rbspb,'Color',colors(j,:),'LineStyle','-.','LineWidth',2)
    elseif strcmp(separate_or_combine,'comb')
        data = [[epoch_rbspa;epoch_rbspb],[TRBEC_rbspa;TRBEC_rbspb]];
        data = sortrows(data,1);
        TRBEC = data(:,2);
        epoch = data(:,1);
        
        TRBEC = smooth(epoch,TRBEC,5);
        
        semilogy(epoch,TRBEC,'Color',colors(j,:),'LineWidth',2)
        hold on
        
    end
end
output_size = 2400*[3.2 1];
resolution = 500;
set(gcf,'Units','inches','Position',[0 0 output_size/resolution]);
end


%given filesnames of the differential content and integration bounds, computes the TRBEC
function [TRBEC, epoch] = TRBEC_E_PA(directory,filenames,energy_min,energy_max,Lstar_min,Lstar_max)

TRBEC = zeros(length(filenames),1);
epoch = zeros(length(filenames),1);
for i = 1:length(filenames)
    load(strcat(directory,filenames(i)))
    if diff_content == 0
        TRBEC(i) = 0;
    else
        %reduce data based on bounds
        selected_energy_indices = find(energy_interp>=energy_min & energy_interp<=energy_max);
        selected_Lstar_indices = find(Lstar_interp>=Lstar_min & Lstar_interp<=Lstar_max);
        
        %compute the TRBEC using new bounds
        TRBEC(i) = trapz(energy_interp(selected_energy_indices),trapz(Lstar_interp(selected_Lstar_indices),diff_content(selected_Lstar_indices,selected_energy_indices),1),2);
        epoch(i) = epoch_select;
    end
end
end


%given directory, satellite, and start and stop dates, outputs
%a list of files according to start/stop times and the epoch of each file as datenum
%satellite: either 'rbspa' or 'rbspb'
%assumes formatting of ''
function [selected_filenames,selected_file_epoch] = select_data_filenames(data_directory,satellite,start_time,end_time)
%for flux data, select out date from filename
%default formatting: 'DC_E_PA_rbspx_yyyymmdd_HHMMSS.mat'
full_filenames = filelist(strcat(data_directory,'DC_E_PA_',satellite,'*.*'));
%convert from string array (double quotes), to character vector (single quotes), necessary to strip date out of filenames
full_filenames_characters = char(full_filenames);
year = str2num(full_filenames_characters(:,15:18));
month = str2num(full_filenames_characters(:,19:20));
day = str2num(full_filenames_characters(:,21:22));
file_epoch = datenum(year,month,day);

%based on given dates, determines which files should be used and puts them in a list
selected_epoch_indices = find(file_epoch>=start_time & file_epoch<=end_time);%includes all full half orbits that occured within the time period
selected_filenames = full_filenames(selected_epoch_indices);
selected_file_epoch = file_epoch(selected_epoch_indices);

end

%given a directory (which can specify the type of file) returns a list of
%the files names in that folder (extracts the names from the "dir" function)
%example, filelist('data/PSD_rbspa*.*') will load all the files that start with 'PSD_rbspa' in the 'data' folder
function output = filelist(data_directory)

files = dir(data_directory);
number_files = length(files);

file_names = strings(number_files,1);%create empty list
for i = 1:number_files
   file_names(i) = files(i).name;
end
output = file_names;
end