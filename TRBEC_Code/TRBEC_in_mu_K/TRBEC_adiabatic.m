%given filesnames of the differential content and integration bounds, computes the TRBEC
function [TRBEC,epoch] = TRBEC_adiabatic(DC_directory,satellite,selected_half_orbit,mu_min,mu_max,Lstar_min,Lstar_max)

    load(strcat(DC_directory,'DC_mu_K_',satellite,'_',datestr(selected_half_orbit,'yyyymmdd_HHMMSS'),'.mat'))

    %if diff_content is intentionally set to zero (scalar), then set TRBEC to zero
    if diff_content == 0
        TRBEC = 0;
        epoch = epoch_select; 
    else
        %reduce data based on bounds
        selected_mu_indices = find(mu_limits>=mu_min & mu_limits<=mu_max);
        selected_Lstar_indices = find(Lstar_limits>=Lstar_min & Lstar_limits<=Lstar_max);

        %compute the TRBEC using new bounds
        TRBEC = trapz(mu_limits(selected_mu_indices),trapz(Lstar_limits(selected_Lstar_indices),diff_content(selected_mu_indices,selected_Lstar_indices),2),1);
        epoch = epoch_select;
    end
end