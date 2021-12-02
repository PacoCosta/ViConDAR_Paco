% for x=1:size(plane_traj,2)
%     
%     coord_x(:,x)=plane_traj{x}(1,1,1);
% end
% clear uu
% for i=1:length(plane_traj)
%    uu{i}=size(plane_traj{i}) ;
% end
% 
% length(uu{1})

clear coord
for i0=1:size(plane_traj{1},3)
    for i1=1:size(plane_traj{1},2)
        for i2=1:size(plane_traj,2)
            coord{i0,i1}(:,i2)=plane_traj{i2}(:,i1,i0);
        end 
    end
end