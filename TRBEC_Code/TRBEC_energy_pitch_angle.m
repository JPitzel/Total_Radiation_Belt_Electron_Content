%from the differential content, computes the TRBEC
function [TRBEC, epoch] = TRBEC_energy_pitch_angle(DC_directory,satellite,selected_half_orbit,energy_min,energy_max,Lstar_min,Lstar_max)

load(strcat(DC_directory,'DC_E_PA_',satellite,'_',datestr(selected_half_orbit,'yyyymmdd_HHMMSS'),'.mat'));

if diff_content == 0
    TRBEC = 0;
    epoch = epoch_select;
else
%reduce data based on bounds
selected_energy_indices = find(energy_interp>=energy_min & energy_interp<=energy_max);
selected_Lstar_indices = find(Lstar_interp>=Lstar_min & Lstar_interp<=Lstar_max);

%compute the TRBEC using new bounds
TRBEC = trapz(energy_interp(selected_energy_indices),trapz(Lstar_interp(selected_Lstar_indices),diff_content(selected_Lstar_indices,selected_energy_indices),1),2);
epoch = epoch_select;
end
end