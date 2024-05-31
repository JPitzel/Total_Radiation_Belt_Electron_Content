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
