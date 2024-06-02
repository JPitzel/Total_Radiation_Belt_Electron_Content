%Takes in a given time and orbit data and return data for the half orbit that occured during the given time
function [PSD_half_orbit,epoch_half_orbit,mu_half_orbit,K_half_orbit,Lstar_half_orbit,energy_half_orbit,alpha_half_orbit] = half_orbit_select_PSD(PSD,epoch,mu,K,Lstar,energy,alpha,orbtimes,epoch_select)
%from the given time, select out the times that denote the beginning and end of the half orbit
[orbtimes_begin, orbtimes_begin_index]= max((orbtimes<=epoch_select).*orbtimes);
orbtimes_end = orbtimes(orbtimes_begin_index + 1);
%from the beginning and end times select out first and last index within the half orbit
%(NOTE: this assumes increasing values of time which is the case for epoch)
all_epoch_indices = find(epoch>=orbtimes_begin & epoch<orbtimes_end);%include start of orbit, exclude end
orbit_begin_index = all_epoch_indices(1);
orbit_end_index = all_epoch_indices(end);

%remove all data that isn't in the selected orbit based on the indices obtained above
PSD_half_orbit = PSD(:,:,orbit_begin_index:orbit_end_index);       
epoch_half_orbit = epoch(orbit_begin_index:orbit_end_index);     
mu_half_orbit = mu(:,:,orbit_begin_index:orbit_end_index);         
K_half_orbit = K(:,:,orbit_begin_index:orbit_end_index);
Lstar_half_orbit = Lstar(:,:,orbit_begin_index:orbit_end_index);     
energy_half_orbit = energy(:,:,orbit_begin_index:orbit_end_index);
alpha_half_orbit = alpha(:,:,orbit_begin_index:orbit_end_index);
end
