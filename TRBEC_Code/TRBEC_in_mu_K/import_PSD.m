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