%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This matlab script reproduces panels for Figure 1 from:
%   Pitzel et al (2024) "Derivations of the Total Radiation Belt Electron
%   Content", JGR - Space Physics
% showing TRBEC integrated in (mu,K,L*) space in March 2013.
% Script "Example_Fig_1_calculate.m" must be run first to generate the
% required intermediate data files.
% To reproduce the entire figure, run this script 4 times with the
% appropriate bounds for mu by changing N in line:
%   mu_bounds = mu_bounds(N,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
directory = 'DC_mu_K\';
export_filename = 'mu_K_1000_to_2000.jpg';
legend_toggle = 0;
bounds_toggle = 1;
output_toggle = 1;

%start and stop days (inclues half orbits that started on the given dates)
epoch_start = datenum(2013,02,27);
epoch_stop = datenum(2013,03,20);

%for each set of integration bounds, combined data or keep seperate
separate_or_combine = 'sep';%'sep' to keep satellite data separate, 'comb' to combine data together
single_or_multiple = 'sgl'; %'sgl' for individual or 'mult' for multiple
 
%mu bounds [mu_min,mu_max], size must match Lstar_bounds below
mu_bounds = [50, 200;
            200, 500;
            500, 1000;
            1000, 2000;];
mu_bounds = mu_bounds(4,:);

%Lstar bounds [Lstar_min,Lstar_max], size must match mu_bounds above
Lstar_min = 3.5;
Lstar_max = 5;
Lstar_bounds = [Lstar_min,Lstar_max;
                Lstar_min,Lstar_max;
                Lstar_min,Lstar_max;];
            
% Line colors (RGB)
colors = [0.000, 0.447, 0.741;
           0.850, 0.325, 0.098;
           0.929, 0.694, 0.125;	
           0.494, 0.184, 0.556;
           0.466, 0.674, 0.188;
           0.301, 0.745, 0.933;	
           0.635, 0.078, 0.184;];

       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
main(directory,epoch_start,epoch_stop,separate_or_combine,single_or_multiple,colors,mu_bounds,Lstar_bounds)
fontsize = 22;
font_type = 'times';
set(gca,'fontname',font_type)

%axis labels
ylabel('TRBEC')
%ylim([10^23 4*10^25])%for 500 to 1000
ylim([1*10^22 5*10^25])%for 1000 to 2000
xlim([epoch_start epoch_stop]);
datetick('x','dd-mmm','keeplimits')
ax = gca; 
ax.FontSize = fontsize; 
labels = string(ax.XAxis.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
ax.XAxis.TickLabels = labels; % set

if legend_toggle == 1
    legend_plots{1} = plot(nan,'Color',colors(1,:));legend_labels{1} = "MagEIS A";
    legend_plots{2} = plot(nan,'Color',colors(2,:));legend_labels{2} = "MagEIS B";
    dim_legend = [0.448, 0.771, 0.1, 0.1];
    [~, hobj, ~, ~] = legend([legend_plots{:}],legend_labels{:},'location',dim_legend,'FontSize',fontsize);
    hl = findobj(hobj,'type','line');
    set(hl,'LineWidth',2);
end

if bounds_toggle == 1
    %bounds for PA and L*
    dim_bounds = [0.55, 0.6045, 0.3, 0.3];
    %str = {[strcat('mu =',32,num2str(Lstar_bounds(1,1)),32,'to',32,num2str(Lstar_bounds(1,2)))],['K = 0 to K_LC \sqrt(G)'],[strcat('L* =',32,num2str(Lstar_bounds(1,1)),32,'to',32,num2str(Lstar_bounds(1,2)))]};
    str = {strcat('$$\mu = ',num2str(mu_bounds(1,1)),'$$ to $$',num2str(mu_bounds(1,2)),'\,\,(MeV/G)$$'),'$$K = 0$$ to $$K_{LC} \,\, (R_E \sqrt{G})$$','$$L^* = 3.5$$ to $$5$$'};

    annotation('textbox',dim_bounds,'String',str,'FitBoxToText','on','FontSize',fontsize,'BackgroundColor','w','FaceAlpha',1,'FontName',font_type,'interpreter','latex');
end

if output_toggle == 1
    %output
    resolution = 500;
    exportgraphics(gcf,export_filename,'Resolution',resolution)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main(directory,epoch_start,epoch_stop,separate_or_combine,individual_or_multiple,colors,mu_bounds,Lstar_bounds)
    [filenames_rbspa,~] = select_data_filenames(directory,'rbspa',epoch_start,epoch_stop);
    [filenames_rbspb,~] = select_data_filenames(directory,'rbspb',epoch_start,epoch_stop);

    [number_of_bounds,~] = size(mu_bounds);

    for j = 1:number_of_bounds

        [TRBEC_rbspa,epoch_rbspa] = TRBEC_mu_K(directory,filenames_rbspa,mu_bounds(j,1),mu_bounds(j,2),Lstar_bounds(j,1),Lstar_bounds(j,2));
        [TRBEC_rbspb,epoch_rbspb] = TRBEC_mu_K(directory,filenames_rbspb,mu_bounds(j,1),mu_bounds(j,2),Lstar_bounds(j,1),Lstar_bounds(j,2));

        %get rid of TRBEC = 0
        epoch_rbspa(TRBEC_rbspa == 0) = [];epoch_rbspb(TRBEC_rbspb == 0) = [];
        TRBEC_rbspa(TRBEC_rbspa == 0) = [];TRBEC_rbspb(TRBEC_rbspb == 0) = [];
        epoch_rbspa(TRBEC_rbspa > 10^30) = [];epoch_rbspb(TRBEC_rbspb > 10^27) = [];
        TRBEC_rbspa(TRBEC_rbspa > 10^30) = [];TRBEC_rbspb(TRBEC_rbspb > 10^27) = [];


        %remove outliers ('mean' defines outliers as points that are three standard deviations from the mean).
        %[TRBEC_rbspa,rbspa_cleaned_indices] = rmoutliers(TRBEC_rbspa);[TRBEC_rbspb,rbspb_cleaned_indices] = rmoutliers(TRBEC_rbspb);
        %epoch_rbspa = epoch_rbspa(rbspa_cleaned_indices == 0); epoch_rbspb = epoch_rbspb(rbspb_cleaned_indices == 0); 

        if strcmp(individual_or_multiple,'mult') == 1
           figure() 
        end

        if strcmp(separate_or_combine,'sep')
            semilogy(epoch_rbspa,TRBEC_rbspa,'Color',colors(j,:),'LineStyle','-','LineWidth',2)
            hold on
            semilogy(epoch_rbspb,TRBEC_rbspb,'Color',colors(j+1,:),'LineStyle','-','LineWidth',2)
        elseif strcmp(separate_or_combine,'comb')
            data = [[epoch_rbspa;epoch_rbspb],[TRBEC_rbspa;TRBEC_rbspb]];
            data = sortrows(data,1);
            TRBEC = smooth(data(:,2),5);
            epoch = data(:,1);

            semilogy(epoch,TRBEC,'Color',colors(j,:),'LineStyle','-','LineWidth',2)
            hold on

        end
        output_size = 2200*[4 1];
        resolution = 500;
        set(gcf,'Units','inches','Position',[0 0 output_size/resolution]);
    end
end

%given filesnames of the differential content and integration bounds, computes the TRBEC
function [TRBEC, epoch] = TRBEC_mu_K(directory,filenames,mu_min,mu_max,Lstar_min,Lstar_max)
    TRBEC = zeros(length(filenames),1);
    epoch = zeros(length(filenames),1);
    for i = 1:length(filenames)
        load(strcat(directory,filenames(i)))
        if diff_content == 0
            TRBEC(i) = 0;
        else
            %reduce data based on bounds
            selected_mu_indices = find(mu_limits>=mu_min & mu_limits<=mu_max);
            selected_Lstar_indices = find(Lstar_limits>=Lstar_min & Lstar_limits<=Lstar_max);

            %compute the TRBEC using new bounds
            TRBEC(i) = trapz(mu_limits(selected_mu_indices),trapz(Lstar_limits(selected_Lstar_indices),diff_content(selected_mu_indices,selected_Lstar_indices),2),1);
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
    full_filenames = filelist(strcat(data_directory,'DC_mu_K_',satellite,'*.*'));
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


