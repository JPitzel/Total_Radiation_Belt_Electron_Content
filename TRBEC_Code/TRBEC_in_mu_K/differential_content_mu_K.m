%given half orbit and filenames to import data from, outputs the differential content
%as a function of first and third adiabatic invariants (mu,L*)
function differential_content_mu_K(epoch_select,satellite,orbtimes,PSD_directory,flux_directory,DC_directory,filename_PSD,filename_flux)

    %constants
    const.R_E = 6371; %km, Earth radius
    const.E_0 = 0.511; %MeV, electron rest energy
    const.B_0 = 0.311653; %G
    const.mu_E = const.B_0*(const.R_E*10^5)^3; %cm^3 G, magnetic moment of Earth

    %import PSD data from given file list
    PSD = []; epoch = []; mu = []; K = []; Lstar = []; energy = []; alpha = [];
    for PSD_file_idx = 1:length(filename_PSD)
        [PSD_i,epoch_i,mu_i,K_i,Lstar_i,energy_i,alpha_i,datainfo] = import_PSD(strcat(PSD_directory,filename_PSD(PSD_file_idx)));
        PSD = cat(3,PSD,PSD_i);
        epoch = cat(1,epoch,epoch_i);
        mu = cat(3,mu,mu_i);
        K = cat(3,K,K_i);
        Lstar = cat(3,Lstar,Lstar_i);
        energy = cat(3,energy,energy_i);
        alpha = cat(3,alpha,alpha_i);
    end

    %selects out the data for the half orbit that occured during the selected time
    [PSD,epoch,mu,K,Lstar,energy,alpha] = half_orbit_select(PSD,epoch,mu,K,Lstar,energy,alpha,orbtimes,epoch_select);
    [K_length,mu_length,Lstar_length] = size(PSD);

    %import B field data from flux data (its interpolated to be synced with PSD data)
    %(note: if moved from this spot in the code must change alpha = 90 index in the get_B_fields function)
    [Bsat,Beq] = get_B_fields(strcat(flux_directory,filename_flux),epoch,Lstar);

    %add nearest point estimate for the loss cone at small mu
    [PSD,mu,K,Lstar,alpha,energy] = K_loss_cone_small_mu(PSD,mu,K,Lstar,alpha,energy,Bsat,const,1);%using K(alpha = 10)
    [PSD,mu,K,Lstar,alpha,energy] = K_loss_cone_small_mu(PSD,mu,K,Lstar,alpha,energy,Bsat,const,length(alpha(:,1,1))-1);%using K(alpha = 170)

    %add PSD = 0 at the loss cone at large mu
    [PSD,mu,K,Lstar,alpha,energy] = K_loss_cone_large_mu(PSD,mu,K,Lstar,alpha,energy);

    %add in K = 0 particles that cannot be seen by the detector
    [PSD,mu,K,Lstar,alpha,energy] = equatorial_particles(PSD,mu,K,Lstar,alpha,energy,Bsat,Beq);

    %sets K as a placeholder value when the PSD is a placeholder (these points are then ignored in the interpolation)
    K(PSD == 0) = -10^31;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %select bounds and resolution for interpolation (NOTE: K is from 0 to K_LC)
    mu_min = 50;
    mu_max = 2000;
    mu_spacing = 10;
    Lstar_min = 3.5;
    Lstar_max = 5;
    Lstar_spacing = 0.1;
    K_interp_length = 200;
    mu_interp_length = round((mu_max-mu_min)/mu_spacing + 1);
    Lstar_interp_length = length((Lstar_min:Lstar_spacing:Lstar_max));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %interpolates over the region of interest to compute the differential content
    %and saves it to a seperate file
    %use log10(PSD) for interpolation, but first change where PSD = 0 to avoid log(0) issues
    PSD(PSD == 0) = 10^-40;
    logPSD = log10(PSD);

    %finds mean Lstar at each time, used to find closest value later on
    [~,~,Lstar_length] = size(Lstar);
    Lstar_mean = zeros(Lstar_length,1);
    for i = 1:Lstar_length
        Lstar_mean(i) = mean(mean(Lstar(:,:,i)));
    end

    Lstar_limits = linspace(Lstar_min,Lstar_max,Lstar_interp_length);
    mu_limits = logspace(log10(mu_min),log10(mu_max),mu_interp_length);


    diff_content = ones(mu_interp_length,Lstar_interp_length);

    for i = 1:length(Lstar_limits)

        Lstar_slice = Lstar_limits(i);
        [~,nearest_Lstar_index] = min(abs(mean(mean(Lstar-Lstar_slice))));

        K_at_Lstar_slice = K(:,:,nearest_Lstar_index); K_at_Lstar_slice = K_at_Lstar_slice(:);
        mu_at_Lstar_slice = mu(:,:,nearest_Lstar_index); mu_at_Lstar_slice = mu_at_Lstar_slice(:);
        logPSD_at_Lstar_slice = logPSD(:,:,nearest_Lstar_index); logPSD_at_Lstar_slice = logPSD_at_Lstar_slice(:);

        mu_at_Lstar_slice(K_at_Lstar_slice<0) = [];
        logPSD_at_Lstar_slice(K_at_Lstar_slice<0) = [];
        K_at_Lstar_slice(K_at_Lstar_slice<0) = [];
        K_at_Lstar_slice(mu_at_Lstar_slice<0) = [];
        logPSD_at_Lstar_slice( mu_at_Lstar_slice<0) = [];
        mu_at_Lstar_slice(mu_at_Lstar_slice<0) = [];

    % %     plot_toggle = 0;
    % %     if plot_toggle == 1
    % %         figure()
    % %         scatter(K_at_Lstar_slice,mu_at_Lstar_slice,20,10.^logPSD_at_Lstar_slice,'filled')
    % %         title(strcat('interpolated','(L* \approx',num2str(Lstar_slice),')'))
    % %         set(gca,'colorscale','log')
    % %         %xlim([0 2.11020679*Lstar_slice-2.26735643])
    % %         xlim([0 0.5])
    % %         ylim([50 200])
    % %         xlabel('K (R_E \surdG)','Interpreter','tex')
    % %         ylabel('\mu (MeV/G)')
    % %         cb = colorbar;
    % %         cb.Label.String = 'PSD [(c/(cm MeV))^3]';
    % %         caxis([10^-4 10^-2])
    % %     end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %sets values of K interpolation
        K_upper_limit = 2.11020679*Lstar_limits(i)-2.26735643; %K_LC(Lstar) function
        K_limits = logspace(-41,log10(K_upper_limit),K_interp_length);%K range is from 0 (or close it it on a log scale) to K_LC
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %calculate standard deviations and means of the data for use in normalization and reverse process later on
        std_K = std(K_at_Lstar_slice(:));
        std_mu = std(mu_at_Lstar_slice(:));
        mean_K = mean(K_at_Lstar_slice(:));
        mean_mu = mean(mu_at_Lstar_slice(:));
        %normalize K and mu data
        K_at_Lstar_slice = (K_at_Lstar_slice - mean_K)/std_K;
        mu_at_Lstar_slice = (mu_at_Lstar_slice - mean_mu)/std_mu;
        %normalize K and mu interpolation points
        [K_interp,mu_interp] = meshgrid(K_limits,mu_limits);
        K_interp = (K_interp - mean_K)/std_K;
        mu_interp = (mu_interp - mean_mu)/std_mu;

        %interpolate
        F = scatteredInterpolant(K_at_Lstar_slice(:),mu_at_Lstar_slice(:),logPSD_at_Lstar_slice(:),'linear','linear');
        logPSD_interp = F(K_interp,mu_interp);

        %flag if interpolation isn't as expected
        if sum(sum(sum(isnan(logPSD_interp)))) ~= 0
            disp('some PSD components were interpolated as NaN!')
            disp(Lstar_slice)
        end

        %unnormalize interpolated points and recover PSD from log10(PSD)
        K_interp = (K_interp*std_K) + mean_K;
        mu_interp = (mu_interp*std_mu) + mean_mu;
        PSD_interp = 10.^logPSD_interp;


        %another flag if the interpolation isn't as expected
        if max(max(max(isinf(PSD_interp)))) ~= 0
            disp('some PSD components were interpolated as Inf!')
            disp(Lstar_slice)
        end



        %calculates the TRBEC integral argument (neglecting constants)
        TRBEC_int_arg = sqrt(mu_interp).*PSD_interp./(Lstar_slice.^2);

        %compute differential content by integrating over K
        diff_content(:,i) = 8.*sqrt(2).*(pi^2).*(const.E_0^(3/2)).*const.mu_E.*trapz(K_limits,TRBEC_int_arg,2);


        % %     %contour plots of integral argument vs mu and K at varying Lstar
        % %         if plot_toggle == 1
        % %             figure()
        % %             scatter(K_interp(:),mu_interp(:),20,PSD_interp(:),'filled')
        % %             title(strcat('interpolated','(L* \approx',num2str(Lstar_slice),')'))
        % %             set(gca,'colorscale','log')
        % %             %xlim([0 2.11020679*Lstar_slice-2.26735643])
        % %             xlim([0 0.5])
        % %             ylim([50 200])
        % %             xlabel('K (R_E \surdG)','Interpreter','tex')
        % %             ylabel('\mu (MeV/G)')
        % %             cb = colorbar;
        % %             cb.Label.String = 'PSD [(c/(cm MeV))^3]';
        % %             caxis([10^-4 10^-2])
        % %             hold on
        % %             contour(K_interp,mu_interp,PSD_interp,(10^-4:8*10^-4:10^-2),'k')
        % %         end
    end

    %saves the differential content to a seperate file
    mu_limits = mu_limits';Lstar_limits = Lstar_limits';%invert these vectors now, unsure if changing it earlier would affect code
    save(strcat(DC_directory,'DC_mu_K_',satellite,datestr(epoch_select,'_yyyymmdd_HHMMSS'),'.mat'),'diff_content','mu_limits','Lstar_limits','epoch_select')

    %plots after integral over K, integrated integral arguments vs mu and Lstar

    % % if plot_toggle == 1
    % %     figure()
    % %     imagesc(Lstar_limits,mu_limits,diff_content);
    % %     title('after K integral')
    % %     set(gca,'YDir','normal')
    % %     set(gca,'colorscale','log')
    % %     cb = colorbar;
    % %     cb.Label.String = '';
    % %     xlabel('L*')
    % %     ylabel('\mu (MeV/G)')
    % %     cb.Label.String = 'PSD Integrated Over K';
    % % end


end
