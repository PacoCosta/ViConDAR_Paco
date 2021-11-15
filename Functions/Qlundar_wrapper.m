%% Header
%
% Qlundar wrapper. Gets the names fron perm_cell and calculates the TI
%
% F.Costa
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2021



function [TI_mean_lidar,TI_mean_WF]=Qlundar_wrapper(input,perm_cell)

% Obtain the names
n_cases=(length(input.freeInp{1,2})*length(input.freeInp{3,2})*length(input.freeInp{4,2})*length(input.timestep_pat_vec {1})*length(input.timeStep_Measurements{1})) ;
for first=1:n_cases
    cummulative_name=1;
    accum{1,first}=perm_cell.OutNames{first,1};
    for second=1:size(perm_cell.values,1)
        for i = 1:size (perm_cell.values{1,1},2)
            [a{i},v{i}]=(setdiff(perm_cell.values{first,1}{i},perm_cell.values{second,1}{i}));  % identify the seeds
        end
        if ~isempty(a{2}) && isempty(a{1}) && isempty(a{3})&& isempty(a{4})&& isempty(a{5})&& isempty(a{6}) && isempty(a{7})&& isempty(a{8})&& isempty(a{9}) && isempty(a{10}) && isempty(a{11}) && isempty(a{12})
            cummulative_name=1+cummulative_name;
            accum{cummulative_name,first}=perm_cell.OutNames{second,1}; % accumulate in the same column the names with same values but different seed
        end
    end
end

% loading the windfields and calculate the mean for the different seed
for col_names=1:size(accum,2)
    for row_names =1:size(accum,1)
        LoadNam{row_names,col_names}= load([input.LidarOutput_dir accum{row_names,col_names}]);
        TI_DATA_lidar {row_names,col_names}= LoadNam{row_names,col_names}.Output.statistics.U.lidar.TI_mean;
        TI_DATA_fullWF{row_names,col_names}= LoadNam{row_names,col_names}.Output.statistics.U.fullWF.TI_mean;
    end
    TI_mean_lidar {1,col_names} = mean (TI_DATA_lidar {row_names,col_names});
    TI_mean_WF    {1,col_names} = mean (TI_DATA_fullWF{row_names,col_names});
end

g=0;