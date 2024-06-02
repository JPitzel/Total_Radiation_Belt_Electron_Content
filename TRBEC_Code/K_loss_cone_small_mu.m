%adds in points for the loss cone (alpha_LC converted to K_LC) according to steps outlined in overleaf documentation
function [PSD,mu,K,Lstar,alpha,energy] = K_loss_cone_small_mu(PSD,mu,K,Lstar,alpha,energy,Bsat,const,alpha_idx)
[~,mu_length,Lstar_length] = size(PSD);
%step 1
h = 100; %km, choosen altitude at which particles precipitate
alpha_LC = asind((1 + ((1.5).*(h./const.R_E)))./(Lstar(alpha_idx,:,:).^(3/2).*(4-((3./Lstar(alpha_idx,:,:)).*(1+(h./const.R_E)))).^(0.25))); %equation from JF's work
K_LC = 2.11020679*Lstar(alpha_idx,:,:)-2.26735643; %equation from Chris
energy_LC = energy(alpha_idx,:,:);

%step 2
mu_LC = zeros(1,mu_length,Lstar_length);
for i = 1:Lstar_length
    mu_LC(1,:,i) = ((energy_LC(1,:,i).^2+2*const.E_0.*energy_LC(1,:,i)).*(sind(alpha_LC(1,:,i)).^2))./(2*const.E_0.*Bsat(i));
end

Lstar_LC = zeros(1,mu_length,Lstar_length);
for mu_idx = 1:mu_length
    for Lstar_idx = 1:Lstar_length
        if isreal(alpha_LC(1,mu_idx,Lstar_idx))
            Lstar_LC(1,mu_idx,Lstar_idx) = interp1(alpha(:,mu_idx,Lstar_idx),Lstar(:,mu_idx,Lstar_idx),alpha_LC(1,mu_idx,Lstar_idx),'linear','extrap');
        else
            alpha_LC(1,mu_idx,Lstar_idx) = -10^31;
            Lstar_LC(1,mu_idx,Lstar_idx) = -10^31;
        end
    end
end
PSD_LC = zeros(1,mu_length,Lstar_length);
for mu_idx = 1:mu_length
    for Lstar_idx = 1:Lstar_length
        if alpha_LC(1,mu_idx,Lstar_idx) > 0
            PSD_LC(1,mu_idx,Lstar_idx) = PSD(alpha_idx,mu_idx,Lstar_idx);
        else
            PSD_LC(1,mu_idx,Lstar_idx) = 0;
        end
    end
end
%concatenate the loss cone with the RBSP data
PSD = cat(1,PSD,PSD_LC);
mu = cat(1,mu,mu_LC);
K = cat(1,K,K_LC);
Lstar = cat(1,Lstar,Lstar_LC);
alpha = cat(1,alpha,alpha_LC);
energy = cat(1,energy,energy_LC);

end