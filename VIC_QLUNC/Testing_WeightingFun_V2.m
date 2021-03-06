function VFinalTotal_Time=Testing_WeightingFun_V2(input,LOS_points,VFinalTotal_TimeInt2,distanceSlices)
if strcmpi(input.flag_probe_weighting,"mean")
    VFinalTotal_Time = mean(VFinalTotal_TimeInt2,'omitnan');
elseif strcmpi(input.flag_probe_weighting,"gaussian")
    % Introducing Gaussian weights in the performance of probe volume:
    % First we create the weights:
    prob_slices_distance = input.ref_plane_dist+(LOS_points.slicesAv*distanceSlices);
    for ind_points=1:size(VFinalTotal_TimeInt2,2)
        %             for i=1:size(VFinalTotal_TimeInt2{ind_points},2)
        %                 VFinalTotal_TimeInt3=VFinalTotal_TimeInt2{ind_points}(:,i);
        VFinalTotal_TimeInt3=VFinalTotal_TimeInt2(:,ind_points);
        % For a given fwhm (fwhm = 2*Rayleigh length) we calculate the normal distribution:
%         accuracy = 10;
        fwhm     = 2*input.distance_av_space;
        
        %Sigma and fwhm are simply related as follows (by definition):
        sigma          = fwhm/2.355;
        distan = linspace(-5*sigma+input.ref_plane_dist,5*sigma+input.ref_plane_dist,size(VFinalTotal_TimeInt3,1));
        gaussian = (1/(sigma*sqrt(2*pi)))*exp(-0.5*((distan-input.ref_plane_dist)/sigma).^2);
        
        %Find the points above the fwhm
        half_max       = (min(gaussian) + max(gaussian)) / 2;
        probe_distance = distan(gaussian >= half_max); % distances above the fwhm
        
        % Distances at slices queried by the lidar may not match the distances wihin the probe volume. We take the closest points. Other method might be interpolate
        
        gaussian_factor = zeros(1, length(prob_slices_distance));% preallocation
        % Round to the millimeter
        distan         = round(distan,4);
        probe_distance = round(probe_distance,4);
        for ind_dist=1:length(prob_slices_distance) % for all the queried points within the probe volume
            
            query_point = find (distan==prob_slices_distance(ind_dist),1);
            if isempty (query_point) % if the queried point is not in the distances vector created to perform the gaussian weights, take the closest
                [~,ind] =min(abs(probe_distance-prob_slices_distance(ind_dist)));
                gaussian_factor(1,ind_dist) = gaussian (distan==probe_distance(ind)); % gaussian factors for the distances queried by the lidar within the probe length
                %                     gaussian_factor(1,ind_dist) = interp2();
                
            else
                gaussian_factor(1,ind_dist) = gaussian (query_point);
            end
        end
        
        
        % performing weighted mean
        VFinalTotal_TimeInt3_NoNans                  = isnan(VFinalTotal_TimeInt3); %finding nans
        gaussian_factor(VFinalTotal_TimeInt3_NoNans) = nan;
        VFinalTotal_Time    (:,ind_points)       = sum(gaussian_factor'.*VFinalTotal_TimeInt3,'omitnan')/sum(gaussian_factor,'omitnan'); %#ok<AGROW>
        
        
        
        % Using cumulative distribution to perform probabilities
        %             cumulative_probability = cumsum(gaussian)*(distan(2)-distan(1));
        %             for ind_cum=1:length(slices_distance)
        %                 [~,pos_cum]=min(abs(distan-slices_distance(ind_cum)));
        
        %                 index_cum{ind_cum}=[pos_cum,pos_cum+1];
        %             end
        %
        %             for ind_cum1=1:length(index_cum)
        %                 cum_prob(1,ind_cum1) =  (distan(index_cum{ind_cum1}(2))*cumulative_probability(index_cum{ind_cum1}(2)))-(distan(index_cum{ind_cum1}(1))*cumulative_probability(index_cum{ind_cum1}(1)));
        %             end
        %
        %             % performing weighted mean with cumulative distribution:
        %             VFinalTotal_Time2(:,i) = sum(cum_prob'.*VFinalTotal_TimeInt3,'omitnan');
        
        %             end
        
    end
    
end
end