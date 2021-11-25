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
        interval_of_confidence = input.distance_av_space; 
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
    for ind_points = 1:size(VFinalTotal_TimeInt2,2)
        VFinalTotal_TimeInt3 = VFinalTotal_TimeInt2(:,ind_points);        
        interval_of_confidence = input.distance_av_space;
        distan     = linspace(-interval_of_confidence,interval_of_confidence,size(VFinalTotal_TimeInt3,1));
        
        %Remove Nans
        VFinalTotal_TimeInt3_NoNans = isnan(VFinalTotal_TimeInt3); %finding nans
        
        % Pulsed lidar Weighting function:
        c = 2.99792458e8; % speed of light        
        Pulsed_WeightFun                              = ((1/(input.tau_meas*c))*(erf((4*sqrt(log(2))*(distan)/((c*input.tau)))+(sqrt(log(2)))*input.tau_meas/input.tau)-erf((4*sqrt(log(2))*(distan)/((c*input.tau)))-(sqrt(log(2)))*input.tau_meas/input.tau))); % Taken from literature (see "LEOSPHERE Pulsed Lidar Principles" Cariou J.) 
        Sum_probabilities_RWF                  = sum((distan(2)-distan(1))*Pulsed_WeightFun); %#ok<*NASGU>
        % Performing weighted mean
        Pulsed_WeightFun(VFinalTotal_TimeInt3_NoNans) = nan;
        VFinalTotal_Time (:,ind_points)        = sum(Pulsed_WeightFun'.*VFinalTotal_TimeInt3,'omitnan')/sum(Pulsed_WeightFun,'omitnan'); %#ok<*AGROW>
    end
end
end
