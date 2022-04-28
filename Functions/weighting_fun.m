%% Header
%
% Weighting function:
% Weights the values retrieved by the lidar within the probe volume.
%
% F.Costa
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2021

%--------------------------------------------------------------------------

function VFinalTotal_Time = weighting_fun(input,VFinalTotal_TimeInt2,distanceSlices,distance_X,focus_distances)
distan2 =   (-distance_X/2:distanceSlices:distance_X/2);
distance_2_av=floor(input.distance_av_space/distanceSlices);

[~,center_distan2] = find(distan2== median(distan2));
% Get the truncated distances
distan2_tr = (distan2(center_distan2)-distance_2_av*distanceSlices:distanceSlices:distan2(center_distan2)+distance_2_av*distanceSlices);
for ind_distance=1:size(distan2_tr,2)
    [~,distan2_tr_ind]=ismembertol(distan2_tr(ind_distance),distan2);
    distan2_tr_index(:,ind_distance)=distan2_tr_ind;
end

if strcmpi(input.flag_probe_weighting,"mean")
    VFinalTotal_Time = mean(VFinalTotal_TimeInt2,'omitnan');
elseif strcmpi(input.flag_probe_weighting,"cw")
    
    % For a given fwhm (FWHM = 2*Rayleigh length) we calculate the normal distribution:
    % FWHM  = 2*input.distance_av_space;
    % By definition sigma and fwhm are simply related as follows:
    % sigma = FWHM/(2*sqrt(2*log(2)));
    
    % Introducing CW weights in the calcualtion of probe volume:
    % Create the weights:
    for ind_points = 1:size(VFinalTotal_TimeInt2,2)
        VFinalTotal_TimeInt3 = VFinalTotal_TimeInt2(:,ind_points);
        %         interval_of_confidence = input.distance_av_space;
        %         distan = linspace(-interval_of_confidence,interval_of_confidence,size(VFinalTotal_TimeInt3,1));
        
        
        % Gaussian Weighting function
        
        CW_WeightFun = (1/(pi))*((input.distance_av_space /input.truncation_val)./((input.distance_av_space/input.truncation_val)^2+((distan2).^2)));%(1/(sigma*sqrt(2*pi)))*exp(-0.5*((distan)/sigma).^2);
        CW_WeightFun_tr = CW_WeightFun(distan2_tr_index); %truncate the Weighting function
        
        CW_WeightFun_NoNans  = CW_WeightFun_tr(~isnan(VFinalTotal_TimeInt3));
        distan2_tr_NoNans    = distan2_tr (~isnan(VFinalTotal_TimeInt3));
        
        VFinalTotal_Time_NoNans_tr = VFinalTotal_TimeInt3(~isnan(VFinalTotal_TimeInt3)); %finding nans
        %         ind_points
        %        [~,nan_ind] = find(isnan(VFinalTotal_TimeInt3'));
        %         FWHM_CW = fwhm(distan,CW_WeightFun);
        
        % Find the points that correspond to the no-nans values in VFinalTotal_TimeInt3_NoNans
        % CW_WeightFun_NoNans  = CW_WeightFun(~isnan(VFinalTotal_TimeInt3));
        %         CW_WeightFun_tr(nan_ind)=[];
        %         CW_WeightFun_NoNans  = CW_WeightFun_tr(~isnan(VFinalTotal_Time_NoNans_tr));
        %         CW_WeightFun_NoNans  = CW_WeightFun_tr;
        %         distan2_tr(nan_ind)=[];
        %         distan_NoNans_tr=distan2_tr(~isnan(VFinalTotal_Time_NoNans_tr));
        %         distan_NoNans_tr=distan2_tr;
        
        % Check that sum of probabilities i
        %         Sum_probabilities_CW = sum((distan(2)-distan(1))*CW_WeightFun); %#ok<*NASGU>
        
        VFinalTotal_Time(:,ind_points) = V_spec(input,VFinalTotal_Time_NoNans_tr,CW_WeightFun_NoNans,distan2_tr_NoNans); %#ok<*AGROW>
        %         VFinalTotal_Time(:,ind_points) = sum(CW_WeightFun'.*VFinalTotal_TimeInt3,'omitnan')/sum(CW_WeightFun,'omitnan'); %#ok<*AGROW>
    end
    
    %
    %     % For a given fwhm (FWHM = 2*Rayleigh length) we calculate the normal distribution:
    %     % FWHM  = 2*input.distance_av_space;
    %     % By definition sigma and fwhm are simply related as follows:
    %     % sigma = FWHM/(2*sqrt(2*log(2)));
    %
    %     % Introducing CW weights in the calcualtion of probe volume:
    %     % First we create the weights:
    %     for ind_points = 1:size(VFinalTotal_TimeInt2,2)
    %         VFinalTotal_TimeInt3 = VFinalTotal_TimeInt2(:,ind_points);
    %         interval_of_confidence = input.distance_av_space;
    %         distan     = linspace(-interval_of_confidence,interval_of_confidence,size(VFinalTotal_TimeInt3,1));
    %         %Remove Nans
    %         VFinalTotal_TimeInt3_NoNans = VFinalTotal_TimeInt3(~isnan(VFinalTotal_TimeInt3)); %finding nans
    %
    %         % Gaussian Weighting function
    %
    %         CW_WeightFun = (1/pi)*((input.distance_av_space/input.truncation_val)./((input.distance_av_space/input.truncation_val)^2+((distan).^2)));%cc=(1/(sigma*sqrt(2*pi)))*exp(-0.5*((distan)/sigma).^2);
    % %         FWHM_CW = fwhm(distan,CW_WeightFun);
    % % FWHM_CW2 = fwhm(distan,CW_WeightFun2);
    %         CW_WeightFun_NoNans  = CW_WeightFun(~isnan(VFinalTotal_TimeInt3));
    %         % Check that sum of probabilities is ~ 0.76 (FWHM)
    %         Sum_probabilities_CW = sum((distan(2)-distan(1))*CW_WeightFun); %#ok<*NASGU>
    %         % Performing weighted mean
    % %         CW_WeightFun(VFinalTotal_TimeInt3_NoNans) = nan;
    %         VFinalTotal_Time(:,ind_points) = V_spec(input,VFinalTotal_TimeInt3_NoNans,CW_WeightFun_NoNans,distan); %#ok<*AGROW>
    % %         VFinalTotal_Time(:,ind_points) = sum(CW_WeightFun'.*VFinalTotal_TimeInt3,'omitnan')/sum(CW_WeightFun,'omitnan'); %#ok<*AGROW>
    %     end
elseif strcmpi(input.flag_probe_weighting,"pulsed")
    for ind_points = 1:size(VFinalTotal_TimeInt2,2)
        VFinalTotal_TimeInt3   = VFinalTotal_TimeInt2(:,ind_points);
        interval_of_confidence = input.distance_av_space;
        distan                 = linspace(-interval_of_confidence,interval_of_confidence,size(VFinalTotal_TimeInt3,1));
        %Remove Nans
        VFinalTotal_TimeInt3_NoNans = VFinalTotal_TimeInt3(~isnan(VFinalTotal_TimeInt3)); %finding nans
        % Pulsed lidar Weighting function:
        c = 2.99792458e8; % speed of light
        Pulsed_WeightFun         = ((1/(input.tau_meas*c))*(erf((4*sqrt(log(2))*(distan)/((c*input.tau)))+(sqrt(log(2)))*input.tau_meas/input.tau)-erf((4*sqrt(log(2))*(distan)/((c*input.tau)))-(sqrt(log(2)))*input.tau_meas/input.tau))); % Taken from literature (see "LEOSPHERE Pulsed Lidar Principles" Cariou J.)
        Sum_probabilities_pulsed = sum((distan(2)-distan(1))*Pulsed_WeightFun); %#ok<*NASGU>
        Pulsed_WeightFun_NoNans  = Pulsed_WeightFun(~isnan(VFinalTotal_TimeInt3));
        %         cc=cumsum(Pulsed_WeightFun)
        %         VFinalTotal_Time (:,ind_points)        = sum(Pulsed_WeightFun'.*VFinalTotal_TimeInt3,'omitnan')/sum(Pulsed_WeightFun,'omitnan'); %#ok<*AGROW>
        %         FWHM_Pulsed= fwhm(distan,Pulsed_WeightFun);
        % Velocity spectra and peak detection methods
        VFinalTotal_Time (:,ind_points) = V_spec(input,VFinalTotal_TimeInt3_NoNans,Pulsed_WeightFun_NoNans,distan);
        ind_points;
    end
end
end
