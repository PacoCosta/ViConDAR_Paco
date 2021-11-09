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
        interval_of_confidence = 12*input.distance_av_space;%4*sigma;% % Taken from literature (see "Comparison of methods to derive radial wind speed from a continuous-wave coherent lidar Doppler spectrum"  Held D. and Mann J. - 2018)
%         if 1==1
%             % if want to weight only the points within the probe volume. Then ~100% of the data (4sigma) is within the probe volume
%             interval_of_confidence = 12*sigma;
%         else
%             % Here we capture the selected points within the probe volume
%             % (input.points_av_slice) taking into account that the
%             % weighting function applies for all the space along the beam
%             interval_of_confidence = fwhm/2;
%         end
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
%         if 1==1
%             % if want to weight only the points within the probe volume. Then ~100% of the data (4sigma) is within the probe volume
%             interval_of_confidence = 4*sigma;
%         else
%             % Here we capture the selected points within the probe volume
%             % (input.points_av_slice) taking into account that the
%             % weighting function applies for all the space along the beam
%             interval_of_confidence = fwhm/2;
%         end
        interval_of_confidence = 12*input.distance_av_space;%4*sigma; % 4sigma implies almost 100% of the data
        distan     = linspace(-interval_of_confidence,interval_of_confidence,size(VFinalTotal_TimeInt3,1));
        
        %Remove Nans
        VFinalTotal_TimeInt3_NoNans = isnan(VFinalTotal_TimeInt3); %finding nans
        
        % Pulsed lidar Weighting function:
        % Parameters to define the pulse shape (These should be placed elsewhere)
%         tau_meas = 665e-9; % Time of measurement
%         tau      = 275e-9; % Define the pulse
        c        = 2.99792458e8;
        
        Pulsed_WF                              = ((1/(input.tau_meas*c))*(erf((4*sqrt(log(2))*(distan)/((c*input.tau)))+(sqrt(log(2)))*input.tau_meas/input.tau)-erf((4*sqrt(log(2))*(distan)/((c*input.tau)))-(sqrt(log(2)))*input.tau_meas/input.tau))); % Taken from literature (see "LEOSPHERE Pulsed Lidar Principles" Cariou J.) 
        Sum_probabilities_RWF                  = sum((distan(2)-distan(1))*Pulsed_WF); %#ok<*NASGU>
        % Performing weighted mean
        Pulsed_WF(VFinalTotal_TimeInt3_NoNans) = nan;
        VFinalTotal_Time (:,ind_points)        = sum(Pulsed_WF'.*VFinalTotal_TimeInt3,'omitnan')/sum(Pulsed_WF,'omitnan'); %#ok<*AGROW>
    end
end
end
