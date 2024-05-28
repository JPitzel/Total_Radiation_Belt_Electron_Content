%adds in PSD = 0 at the loss cone for large mu
function [PSD,mu,K,Lstar,alpha,energy] = K_loss_cone_large_mu(PSD,mu,K,Lstar,alpha,energy)
[~,mu_length,Lstar_length] = size(PSD);
%add in K_LC for high values of mu (make interpolatation easier and hence the integral from 0 to K_LC
PSD = cat(1,PSD,10^-40*ones(1,mu_length,Lstar_length));%use small number for PSD and not zero because PSD = 0 is fill value of data
mu = cat(1,mu,3000*ones(1,mu_length,Lstar_length));
K = cat(1,K,K(end,:,:));
Lstar = cat(1,Lstar,Lstar(end,:,:));
alpha = cat(1,alpha,zeros(1,mu_length,Lstar_length));%values are meaningless, but need to something so dimensions are consistent between every variable
energy = cat(1,energy,zeros(1,mu_length,Lstar_length));%these are also meaningless
end
