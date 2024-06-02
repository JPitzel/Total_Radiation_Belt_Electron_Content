%takes in data and selects out data for a given half orbit
function [flux_half_orbit,epoch_half_orbit,Lstar_half_orbit,B_sat_half_orbit,B_eq_half_orbit] = half_orbit_select(flux,epoch,Lstar,B_sat,B_eq,orbtimes,epoch_select)
%from the given time, select out the times at the beginning and end of the half orbit
[orbtimes_begin, orbtimes_begin_index]= max((orbtimes<=epoch_select).*orbtimes);
orbtimes_end = orbtimes(orbtimes_begin_index + 1);
%from the beginning and end times select out first and last index within the half orbit
%(NOTE: this assumes increasing values of time which is the case for epoch)
all_epoch_indices = find(epoch>=orbtimes_begin & epoch<orbtimes_end);%include start of orbit, exclude end (since it's the beginning of the next)
orbit_begin_index = all_epoch_indices(1);
orbit_end_index = all_epoch_indices(end);


%remove all data that isn't in the selected orbit based on the indices obtained above
flux_half_orbit = flux(:,:,orbit_begin_index:orbit_end_index);
epoch_half_orbit = epoch(orbit_begin_index:orbit_end_index);
Lstar_half_orbit = Lstar(orbit_begin_index:orbit_end_index);
B_sat_half_orbit = B_sat(orbit_begin_index:orbit_end_index);
B_eq_half_orbit = B_eq(orbit_begin_index:orbit_end_index);

%remove all data that doesn't have an associated L* value
flux_half_orbit(:,:,Lstar_half_orbit<0) = [];
epoch_half_orbit(Lstar_half_orbit<0) = [];
B_sat_half_orbit(Lstar_half_orbit<0) = [];
B_eq_half_orbit(Lstar_half_orbit<0) = [];
Lstar_half_orbit(Lstar_half_orbit<0) = [];
end
