%adds in PSD = 0 and ghost points for interpolation in mu/K. This ensures interpolation functions properly.
function [PSD,mu,K,Lstar,alpha,energy] = K_mu_ghost_points(PSD,mu,K,Lstar,alpha,energy)
[~,mu_length,Lstar_length] = size(PSD);


%add in K_LC for high values of mu (make interpolatation easier and hence the integral from 0 to K_LC
PSD = cat(1,PSD,10^-40*ones(1,mu_length,Lstar_length));%use small number for PSD and not zero because PSD = 0 is fill value of data
mu = cat(1,mu,11000*ones(1,mu_length,Lstar_length));
K = cat(1,K,K(end,:,:));
Lstar = cat(1,Lstar,Lstar(end,:,:));
alpha = cat(1,alpha,zeros(1,mu_length,Lstar_length));%These are to match shapes
energy = cat(1,energy,zeros(1,mu_length,Lstar_length));%these are to match shapes

PSD = cat(1,PSD,10^-40*ones(1,mu_length,Lstar_length));%use small number for PSD and not zero because PSD = 0 is fill value of data
mu = cat(1,mu,zeros(1,mu_length,Lstar_length));
K = cat(1,K,10*ones(1,mu_length,Lstar_length));
Lstar = cat(1,Lstar,Lstar(end,:,:));
alpha = cat(1,alpha,zeros(1,mu_length,Lstar_length));%These are to match shapes
energy = cat(1,energy,zeros(1,mu_length,Lstar_length));%these are to match shapes

PSD = cat(1,10^-40*ones(1,mu_length,Lstar_length),PSD);
mu = cat(1, 11000*ones(1,mu_length,Lstar_length), mu);
K  = cat(1,zeros(1,mu_length,Lstar_length), K);
Lstar = cat(1,Lstar(1,:,:),Lstar);
alpha = cat(1,zeros(1,mu_length,Lstar_length),alpha);%These are to match shapes
energy = cat(1,zeros(1,mu_length,Lstar_length), energy);%These are to match shapes
end
