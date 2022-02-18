%% Header
%
% Interpolate to get all points (even when we are not in the grid) taking
% into account LOS and position of LiDAR.
%
% V.Pettas/F.Costa
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2019
%--------------------------------------------------------------------------

function [VFinalTotal,VFinalTotal_Time,Y1,Z1] = interpolationFun(input,component,LOS_points,gridy,gridz,distanceSlices,slicesDistance,focus_distances,points_probeX)
if strcmpi(input.interpolation_slices,'interpolate')
    %     for i1 = 1:size(LOS_points.slices,1) % loop over the points of pattern
    %         Y1 = LOS_points.Coor{i1}(1,:);
    %         Z1 = LOS_points.Coor{i1}(2,:);
    %         for ind_slice=1:size(LOS_points.slices,2)
    %             for iTSlice = 1:length(LOS_points.slicesAv) % For measured slices in the pattern
    %                 indLoopT  =  (LOS_points.slices(i1,ind_slice)-1)*distanceSlices+(focus_distances(iTSlice)-input.ref_plane_dist);% [m] Distances where the measurements are focused  along the probe length
    %                 indLoopT2 = indLoopT;
    %                 indNEg = find(indLoopT<=0); % find negative, zeros or Nans
    %                 indNEg =  [indNEg find(isnan(indLoopT))]; % find negative or Nans
    %                 indNEg = [indNEg find(indLoopT>size(component,2)*distanceSlices)]; % findd points outside of the grid (we assume squared grid)
    %                 indLoopT2(indNEg) = [];
    %                 maxLoopInd = round(length(indLoopT)/2); % find the middle of the total slices
    %                 NansStart  = length(find(indNEg<maxLoopInd));
    %                 NansEnd    = length(indNEg)-NansStart;
    %                 indLoopT2=[nan(1,NansStart) indLoopT2 nan(1,NansEnd) ];                     %#ok<*AGROW>
    %                 % Query points for interpolation:
    %                 xq1 = LOS_points.Coor{i1}(1,iTSlice);
    %                 xq2 = LOS_points.Coor{i1}(2,iTSlice);
    %                 xq3 = indLoopT2;
    %                 % Interpolation
    %                 VFinalTotal_TimeInt2{i1}(iTSlice,ind_slice)=interpn(gridz,slicesDistance,gridy,component,xq2,xq3,xq1);
    %             end
    %         end
    %         if length(LOS_points.slicesAv) ~= 1
    %             VFinalTotal_TimeInt3=VFinalTotal_TimeInt2{i1};
    %             VFinalTotal_Time{i1} = weighting_fun(input,LOS_points,VFinalTotal_TimeInt3,distanceSlices);
    %         else
    %             VFinalTotal_Time{i1} = VFinalTotal_TimeInt2;
    %         end
    %     end
    %     % For the complete WF, insead of taking the closest, interpolate along time (x axis):
    %     for ind_s=1:size(LOS_points.Coor,2)
    %         Y2=LOS_points.Coor{ind_s}(1,:);
    %         Z2=LOS_points.Coor{ind_s}(2,:);
    %         Mid_Y2=Y2(floor(length(Y2)/2)+1); % find the point in the middle of the vector
    %         Mid_Z2=Z2(floor(length(Z2)/2)+1);
    %         VFinalTotal{ind_s} = interpn(gridz,slicesDistance,gridy,component,Mid_Z2,slicesDistance,Mid_Y2);
    %     end
    for i1 = 1:size(LOS_points.slices,1) % loop over the points of pattern
        X1 = LOS_points.Coor{i1}(1,:);
        Y1 = LOS_points.Coor{i1}(2,:);
        Z1 = LOS_points.Coor{i1}(3,:);
        Y1 = round(Y1,5); % Round because if not, the interpolation gives NaN's or incorrect results
        Z1 = round(Z1,5);
        X1 = round(X1,5);
        [~,zeros_X_ind]= find(X1<0);
        X1(zeros_X_ind)=[];
        Y1(zeros_X_ind)=[];
        Z1(zeros_X_ind)=[];
        busquedaX = find(ismember(slicesDistance,X1) ); % look for coincidences in X component
        busquedaY = find(ismember(gridy,Y1) ); % look for coincidences in Y component
        busquedaZ = find(ismember(gridz,Z1)); %#ok<*EFIND> % look for coincidences in  Z component
        
        if isempty (busquedaX) || length(busquedaX) <(length(LOS_points.slicesAv)) %check if all measured points are grid points
            for i = 1:length(Y1)    % find the closest points in the grid and use these
                
                DifX = slicesDistance-X1(i);
                [~,indFX] = min(abs(DifX));
                PointFm{i1}(1,i) = slicesDistance(indFX);
                PoinInd{i1}(1,i) = indFX;
                DifX(indFX)=nan;
                [~,indFX] = min(abs(DifX));
                PointFm{i1}(2,i) = slicesDistance(indFX);
                PoinInd{i1}(2,i) = indFX;
                
            end
        else
            PointFm{i1}(1,:) = gridy(busquedaX);
            PoinInd{i1}(1,:) = busquedaX;
            PointFm{i1}(2,:) = gridy(busquedaX);
            PoinInd{i1}(2,:) = busquedaX;
        end
        
        if isempty (busquedaY) || length(busquedaY) <(length(LOS_points.slicesAv)) %check if all measured points are grid points
            for i = 1:length(Y1)    % find the closest points in the grid and use these
                
                DifY = gridy-Y1(i);
                [~,indFY] = min(abs(DifY));
                PointFm{i1}(3,i) = gridy(indFY);
                PoinInd{i1}(3,i) = indFY;
                DifY(indFY)=nan;
                [~,indFY] = min(abs(DifY));
                PointFm{i1}(4,i) = gridy(indFY);
                PoinInd{i1}(4,i) = indFY;
            end
        else
            PointFm{i1}(3,:) = gridy(busquedaY);
            PoinInd{i1}(3,:) = busquedaY;
            PointFm{i1}(4,:) = gridy(busquedaY);
            PoinInd{i1}(4,:) = busquedaY;
        end
        if isempty (busquedaZ) || length(busquedaZ) <(length(LOS_points.slicesAv)) %check if all measured points are grid points
            for i = 1:length(Y1)
                DifZ = gridz-Z1(i);
                [~,indF] = min(abs(DifZ));
                PointFm{i1}(5,i) = gridz(indF);
                PoinInd{i1}(5,i) = indF;
                DifZ(indF)=nan;
                [~,indF] = min(abs(DifZ));
                PointFm{i1}(6,i) = gridz(indF);
                PoinInd{i1}(6,i) = indF;
            end
        else
            PointFm{i1}(5,:) = gridz(busquedaZ); % points in the grid in meters matching our points (nearest, not exactly the value of the trajectory!!!)
            PoinInd{i1}(5,:) = busquedaZ; %indices of point matching in the grid
            PointFm{i1}(6,:) = gridz(busquedaZ); % points in the grid in meters matching our points (nearest, not exactly the value of the trajectory!!!)
            PoinInd{i1}(6,:) = busquedaZ; %indices of point matching in the grid
            
        end
        for ind_col=1:size(PointFm{i1},2)
            points2interp_ind{ind_col}  = combvec(PoinInd{i1}(1:2,ind_col)',PoinInd{i1}(3:4,ind_col)',PoinInd{i1}(5:6,ind_col)');
            points2interp_val {ind_col} = combvec(PointFm{i1}(1:2,ind_col)',PointFm{i1}(3:4,ind_col)',PointFm{i1}(5:6,ind_col)');
            
            for i_iint=1:size(points2interp_val{ind_col},2)
                v_points2interp{ind_col}(:,i_iint) = component(points2interp_ind{ind_col}(3,i_iint),points2interp_ind{ind_col}(1,i_iint),points2interp_ind{ind_col}(2,i_iint));
            end
            resh_vel_val{ind_col} = reshape(v_points2interp{ind_col},2,2,2);
            int_gridx{ind_col}    = unique(points2interp_val{ind_col}(1,:));
            int_gridy{ind_col}    = unique(points2interp_val{ind_col}(2,:));
            int_gridz{ind_col}    = unique(points2interp_val{ind_col}(3,:));
            
            VFinalTotal_Time{i1}(:,ind_col)=interpn(int_gridz{ind_col},int_gridx{ind_col},int_gridy{ind_col},resh_vel_val{ind_col},Z1(ind_col),X1(ind_col),Y1(ind_col));
        end

        % For the complete WF, insead of taking the closest, interpolate along time (x axis):
        for ind_s=1:size(LOS_points.Coor,2)
            Y2=LOS_points.Coor{ind_s}(1,:);
            Z2=LOS_points.Coor{ind_s}(2,:);
            Mid_Y2=Y2(floor(length(Y2)/2)+1); % find the point in the middle of the vector
            Mid_Z2=Z2(floor(length(Z2)/2)+1);
            VFinalTotal{ind_s} = interpn(gridz,slicesDistance,gridy,component,Mid_Z2,slicesDistance,Mid_Y2);
        end
        
        
        int_gridx2=[0,1];
        int_gridy2=[0,1];
        int_gridz2=[0,1];
        C=reshape([0,1,2,3,4,5,6,7],2,2,2);
        
        V_final2=interpn(int_gridx2,int_gridy2,int_gridz2,C,1,1,1);
        
    end
elseif strcmpi(input.interpolation_slices,'none') %if you don't interpolate get the closest point
    for i1 = 1:size(LOS_points.slices,1) % loop over the points of pattern
        Y1 = LOS_points.Coor{i1}(2,:);
        Z1 = LOS_points.Coor{i1}(3,:);
        Y1 = round(Y1,5); % Round because if not, the interpolation gives NaN's or incorrect results
        Z1 = round(Z1,5);
        busquedaY = find(ismember(gridy,Y1) ); % look for coincidences in Y component
        busquedaZ = find(ismember(gridz,Z1)); %#ok<*EFIND> % look for coincidences in  Z component
        if isempty (busquedaY) || length(busquedaY) <(length(LOS_points.slicesAv)) %check if all measured points are grid points
            for i = 1:length(Y1)    % find the closest points in the grid and use these
                DifY = gridy-Y1(i);
                [~,indF] = min(abs(DifY));
                PointFm{i1}(1,i) = gridy(indF);
                PoinInd{i1}(1,i) = indF;
            end
        else
            PointFm{i1}(1,:) = gridy(busquedaY);
            PoinInd{i1}(1,:) = busquedaY;
        end
        if isempty (busquedaZ) || length(busquedaZ) <(length(LOS_points.slicesAv)) %check if all measured points are grid points
            for i = 1:length(Y1)
                DifZ = gridz-Z1(i);
                [~,indF] = min(abs(DifZ));
                PointFm{i1}(2,i) = gridz(indF);
                PoinInd{i1}(2,i) = indF;
            end
        else
            PointFm{i1}(2,:) = gridz(busquedaZ); % points in the grid in meters matching our points (nearest, not exactly the value of the trajectory!!!)
            PoinInd{i1}(2,:) = busquedaZ; %indices of point matching in the grid
        end
        % Now take the values at this point
        if length(LOS_points.slicesAv)~=1
            VFinalTotal{i1} = squeeze(component(PoinInd{i1}(2,floor(length(LOS_points.slicesAv)/2+1)),:,PoinInd{i1}(1,floor(length(LOS_points.slicesAv)/2+1)))); % the floor is done to always take the middle point of consecuteive measurement ranges
            VFinalTotal{i1} = round(VFinalTotal{i1},5);
        else
            VFinalTotal{i1} = squeeze(component(PoinInd{i1}(2,1),:,PoinInd{i1}(1,1))); % the floor is done to always take the middle point of consecuteive measurement ranges
            VFinalTotal{i1} = round(VFinalTotal{i1},5);
        end
        for iTSlice = 1:length(LOS_points.slicesAv)% loop over the ranges of each point of the pattern and average them
            indLoopT =  LOS_points.slices(i1,:)+LOS_points.slicesAv(iTSlice);
            indLoopT2 = indLoopT;
            indNEg = find(indLoopT<=0); % find negative, zeros or Nans
            indNEg =  [indNEg find(isnan(indLoopT))]; % find negative or Nans
            indNEg = [indNEg find(indLoopT>size(component,2))]; % findd points outside of the grid
            indLoopT2(indNEg) = [];
            VFinalTotal_TimeInt{iTSlice} = squeeze(component(PoinInd{i1}(2,iTSlice), indLoopT2 ,PoinInd{i1}(1,iTSlice)));    %
            % find how many nan points to add in the beginning or in the
            % end to keep length consistent
            maxLoopInd = round(length(indLoopT)/2); % find the middle of the total slices
            NansStart  = length(find(indNEg<maxLoopInd));
            NansEnd    = length(indNEg)-NansStart;
            VFinalTotal_TimeInt2(iTSlice,:) = [nan(1,NansStart) VFinalTotal_TimeInt{iTSlice} nan(1,NansEnd) ];
        end
        if length(LOS_points.slicesAv) ~= 1
            %             VFinalTotal_Time{i1} = weighting_fun(input,LOS_points,VFinalTotal_TimeInt2,distanceSlices);
            VFinalTotal_Time{i1} = weighting_fun(input,VFinalTotal_TimeInt2);
        else
            VFinalTotal_Time{i1} = VFinalTotal_TimeInt2;
        end
        clear VFinalTotal_TimeInt
    end
end
