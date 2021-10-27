%% Header
%
% Weighting function:
% Weights the values retrieved by the lidar within the probe volume.
%
% V.Pettas/F.Costa
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2021

%--------------------------------------------------------------------------

function VFinalTotal_Time=Testing_WeightingFun_V2(input,LOS_points,VFinalTotal_TimeInt2,distanceSlices)
if strcmpi(input.flag_probe_weighting,"mean")
    VFinalTotal_Time = mean(VFinalTotal_TimeInt2,'omitnan');
elseif strcmpi(input.flag_probe_weighting,"gaussian")    
    % Introducing Gaussian weights in the performance of probe volume:
    % First we create the weights:    
    for ind_points = 1:size(VFinalTotal_TimeInt2,2)
        VFinalTotal_TimeInt3 = VFinalTotal_TimeInt2(:,ind_points);        
        % For a given fwhm (fwhm = 2*Rayleigh length) we calculate the normal distribution:
        fwhm  = 2*input.distance_av_space;        
        % By definition sigma and fwhm are simply related as follows:
        sigma = fwhm/(2*sqrt(2*log(2)));
        if 1==1 % in the future want to implement the possibility of assessing CW lidars. The interval we apply the weighting function changes depending on the type of lidar
            interval_of_confidence = 4*sigma; % limits are chosen to ensure a high interval of confidence (4 sigma). This limit is for pulsed lidar, since we catch all the points within the pulse
        else
            interval_of_confidence = input.distance_av_space; % limits for CW, since we all catch the points within the probe length
        end
        distan                 = linspace(-interval_of_confidence,interval_of_confidence,size(VFinalTotal_TimeInt3,1));
        gaussian_w             = (1/(sigma*sqrt(2*pi)))*exp(-0.5*((distan)/sigma).^2);        
        % Check that sum of probabilities is ~ 1.0
        Sum_probabilities=sum((distan(2)-distan(1))*gaussian_w); %#ok<*NASGU>       
        % Performing weighted mean
        VFinalTotal_TimeInt3_NoNans             = isnan(VFinalTotal_TimeInt3); %finding nans
        gaussian_w(VFinalTotal_TimeInt3_NoNans) = nan;
        VFinalTotal_Time    (:,ind_points)      = sum(gaussian_w'.*VFinalTotal_TimeInt3,'omitnan')/sum(gaussian_w,'omitnan'); %#ok<*AGROW>
    end
end
end
