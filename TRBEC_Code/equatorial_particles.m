%add in K = 0 particles that cannot be seen by the detector
%(NOTE: small K corresponds to pitch angles close to 90 so index 17 is alpha = 90 which is the smallest K)
function [PSD,mu,K,Lstar,alpha,energy] = equatorial_particles(PSD,mu,K,Lstar,alpha,energy,Bsat,Beq)
[~,mu_length,Lstar_length] = size(PSD);
K = cat(1,zeros(1,mu_length,Lstar_length),K);
mu = cat(1,mu(17,:,:),mu);
Lstar = cat(1,Lstar(17,:,:),Lstar);
PSD = cat(1,PSD(17,:,:),PSD);
%PSD = cat(1,zeros(1,mu_length,Lstar_length),PSD);%sets PSD = 0 for testing, switch with line above
alpha = cat(1,90*ones(1,mu_length,Lstar_length),alpha);%K = 0 corresponds to alpha = 90;
energy = cat(1,energy(17,:,:),energy);%not sure, but energies should be close to the energies at 90

%corrects all mu values at K = 0 from nearest point estimate to better
%approximation using curves of constant energy
for i = 1:Lstar_length
    mu(1,:,i) = Bsat(i)./Beq(i).*mu(1,:,i);
end
end