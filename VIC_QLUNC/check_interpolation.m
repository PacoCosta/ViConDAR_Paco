


for indexx=1:size(Y1,2)
    xq1 = Y1(1,indexx);
    xq2 = Z1(1,indexx);
    xq3 = indLoopT2(indexx);               

    % Interpolation
    VFinalTotal_Time4{i1}(1,indexx)=interpn(gridy,input.slicesDistance,gridz,component,xq1,xq3,xq2);  
    
    
end



     % ####### Check the interpolation is working properly#####################               

    %                 busquedaY = find(ismember(gridy,Y1{dim1,dim2}(1,ind_p)) ); % look for coincidences in Y component
    %                 busquedaZ = find(ismember(gridz,Z1{dim1,dim2}(1,ind_p))); %#ok<*EFIND> % look for coincidences in  Z component
    %                 if isempty(busquedaY)%|| length(busquedaY) <(length(LOS_points.slicesAv))
    %                     DifY = gridy - Y1{dim1,dim2}(1,ind_p);
    %                     [~,ind_miny] = mink(abs(DifY),2,2);
    %                     ind_miny=sort(ind_miny);
    %                     Point_ind(1,:) = [(ind_miny(1)) (ind_miny(2))];
    %                     Point(1,:) = [gridy(ind_miny(1)) gridy(ind_miny(2))];
    %                 else
    %                     Point_ind(1,:) = [busquedaY,busquedaY];
    %                     Point(1,:)= [gridy(busquedaY) gridy(busquedaY)];
    %                 end
    %                 if isempty(busquedaZ)%|| %length(busquedaZ) <(length(LOS_points.slicesAv))
    %                     DifZ = gridz-Z1{dim1,dim2}(1,ind_p);
    %                     [~,ind_minz] = mink(abs(DifZ),2,2); % we get the two closest points (previous and following ones)
    %                     ind_minz=sort(ind_minz);
    %                     Point_ind(2,:) = [(ind_minz(1)) (ind_minz(2))];
    %                     Point(2,:) = [gridz(ind_minz(1)) gridz(ind_minz(2))];
    %                     
    %                 else
    %                     Point_ind(2,:) = (busquedaZ);
    %                     Point(2,:)= [gridy(busquedaZ) gridy(busquedaZ)];
    %                 end       
    %                 Point_ind(3,:)=input.focus_distances_index{ind_p};
    %                 Point(3,:) = input.focus_distances_new{ind_p};
    %                 point_to_interpolate{ind_p} = (combvec(Point_ind(1,:),Point_ind(3,:),Point_ind(2,:))');
    %                 point_to_interpolate_val{ind_p} = (combvec(Point(1,:),Point(3,:),Point(2,:))');
    %                 % Get velocity values
    %                 for ind_velo_vec=1:size(point_to_interpolate{ind_p},1)
    %                     point_to_interpolate_val{ind_p}(ind_velo_vec,4)= component(point_to_interpolate{ind_p}(ind_velo_vec,1),point_to_interpolate{ind_p}(ind_velo_vec,2),point_to_interpolate{ind_p}(ind_velo_vec,3));
    %                 end
    %                 % interpolate
    %                 point_to_interpolate_val{ind_p}=sortrows(point_to_interpolate_val{ind_p},2);
    %                 x1  = unique(point_to_interpolate_val{ind_p}(:,1)); % Y component
    %                 x2  = flipud(unique(point_to_interpolate_val{ind_p}(:,3))); % Z component
    %                 x3  = unique(point_to_interpolate_val{ind_p}(:,2)); % focus distance (X component)
    %                 Vel_val  = point_to_interpolate_val{ind_p}(:,4); % velocity vector

                    % cases for interpolation: In total there are 8 different
                    % cases to interpolate (have to implement the rest)

    %                 if isempty(busquedaY)==0  && isempty( busquedaZ)==0 && length(x3)~=1 % if Y and Z exist in the grid, but not x interpolate only in the x direction (along the focus distance)
    %                     Vel_mat{1,dim2}(1,ind_p) = interp1(x3, unique(Vel_val,'stable'),xq3,'linear');
    % %                 
    % %                 elseif isempty(busquedaY)==1  && isempty( busquedaZ)==0 
    % %                     Vel_mat{1,dim2}(1,ind_p) = interp2(x1,x3, unique(Vel_val,'stable'),xq1,xq3);
    % %                 
    % %                 elseif isempty( busquedaY)==0  && isempty( busquedaZ)==1
    % %                     Vel_mat{1,dim2}(1,ind_p) = interp2(x2,x3, unique(Vel_val,'stable'),xq2,xq3);
    % %                 
    % %                 elseif isempty(busquedaY)==1  && isempty(busquedaZ)==1 && length(x3)==1
    % %                     Vel_mat{1,dim2}(1,ind_p) = interp2(x1,x2, unique(Vel_val,'stable'),xq1,xq2);                     
    %                 
    %                 elseif isempty( busquedaY)==1  && isempty( busquedaZ)==1 && length(x3)~=1
    %                     V   = reshape(Vel_val,2,2,2);
    %                     Vel_mat{1,dim2}(1,ind_p)= interpn(x1,x2,x3,V,xq1,xq2,xq3,'linear');
    %                 end                                
    % ##########################################################################