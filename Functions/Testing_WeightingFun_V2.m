%% Header
%
% Weighting function:
% Weights the values retrieved by the lidar within the probe volume.
%
% V.Pettas/F.Costa
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2021

%--------------------------------------------------------------------------

function VFinalTotal_Time = Testing_WeightingFun_V2(input,VFinalTotal_TimeInt2)
if strcmpi(input.flag_probe_weighting,"mean")
    VFinalTotal_Time = mean(VFinalTotal_TimeInt2,'omitnan');
elseif strcmpi(input.flag_probe_weighting,"gaussian")
    % For a given fwhm (fwhm = 2*Rayleigh length) we calculate the normal distribution:
    fwhm  = 2*input.distance_av_space;
    % By definition sigma and fwhm are simply related as follows:
    sigma = fwhm/(2*sqrt(2*log(2)));
    
    % Introducing Gaussian weights in the calcualtion of probe volume:
    % First we create the weights:
    for ind_points = 1:size(VFinalTotal_TimeInt2,2)
        VFinalTotal_TimeInt3 = VFinalTotal_TimeInt2(:,ind_points);

        if 1==0
            % if want to weight only the points within the probe volume. Then ~100% of the data (4sigma) is within the probe volume
            interval_of_confidence = 4*sigma;
        else
            % Here we capture the selected points within the probe volume
            % (input.points_av_slice) taking into account that the
            % weighting function applies for all the space along the beam
            interval_of_confidence = fwhm/2;
        end
        distan     = linspace(-interval_of_confidence,interval_of_confidence,size(VFinalTotal_TimeInt3,1));
        
        %Remove Nans
        VFinalTotal_TimeInt3_NoNans = isnan(VFinalTotal_TimeInt3); %finding nans
        
        % Gaussian Weighting function
        gaussian_w = (1/(sigma*sqrt(2*pi)))*exp(-0.5*((distan)/sigma).^2);
        % Check that sum of probabilities is ~ 0.76 (FWHM)
        Sum_probabilities_gaussian = sum((distan(2)-distan(1))*gaussian_w); %#ok<*NASGU>
        % Performing weighted mean
        gaussian_w(VFinalTotal_TimeInt3_NoNans) = nan;
        VFinalTotal_Time    (:,ind_points)      = sum(gaussian_w'.*VFinalTotal_TimeInt3,'omitnan')/sum(gaussian_w,'omitnan'); %#ok<*AGROW>
    end
elseif strcmpi(input.flag_probe_weighting,"pulsed")
    
    % For a given fwhm (fwhm = 2*Rayleigh length) we calculate the normal distribution:
    fwhm  = 2*input.distance_av_space;
    % By definition sigma and fwhm are simply related as follows:
    sigma = fwhm/(2*sqrt(2*log(2)));
    
    for ind_points = 1:size(VFinalTotal_TimeInt2,2)
        VFinalTotal_TimeInt3 = VFinalTotal_TimeInt2(:,ind_points);        
        if 1==0
            % if want to weight only the points within the probe volume. Then ~100% of the data (4sigma) is within the probe volume
            interval_of_confidence = 4*sigma;
        else
            % Here we capture the selected points within the probe volume
            % (input.points_av_slice) taking into account that the
            % weighting function applies for all the space along the beam
            interval_of_confidence = fwhm/2;
        end
        distan     = linspace(-interval_of_confidence,interval_of_confidence,size(VFinalTotal_TimeInt3,1));
        
        %Remove Nans
        VFinalTotal_TimeInt3_NoNans = isnan(VFinalTotal_TimeInt3); %finding nans
        
        % Pulsed lidar Weighting function:
        % Parameters to define the pulse shape (These should be placed elsewhere)
        tau_meas = 265e-9; % Time of measurement
        tau      = 165e-9; % Define the pulse
        c        = 2.99792458e8;
        
        RWF                              = ((1/(tau_meas*c))*(erf((4*sqrt(log(2))*(distan)/((c*tau)))+(sqrt(log(2)))*tau_meas/tau)-erf((4*sqrt(log(2))*(distan)/((c*tau)))-(sqrt(log(2)))*tau_meas/tau)));
        Sum_probabilities_RWF            = sum((distan(2)-distan(1))*RWF); %#ok<*NASGU>
        % Performing weighted mean
        RWF(VFinalTotal_TimeInt3_NoNans) = nan;
        VFinalTotal_Time (:,ind_points)  = sum(RWF'.*VFinalTotal_TimeInt3,'omitnan')/sum(RWF,'omitnan'); %#ok<*AGROW>
    end
end
end
