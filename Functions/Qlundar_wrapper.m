%% Header
%
% Qlundar wrapper. Gets the names fron perm_cell and calculates the mean TI
% from the different seeds
%
% F.Costa
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2021



function TI_Qlundar=Qlundar_wrapper(input,perm_cell)
ii_pat=0;
for ind_pattern=1: length(input.PatternNames)
    
    %     Obtain the names
    n_cases=(length(input.freeInp{1,2})*length(input.freeInp{3,2})*length(input.freeInp{4,2})*length(input.timestep_pat_vec {1})*length(input.timeStep_Measurements{1})) ;
    for cases=1:n_cases
        cummulative_name=1;
        accum{1,cases}=perm_cell.OutNames{cases+ii_pat,1};
        for perm_vals=1:size(perm_cell.values,1)
            for i = 1:size (perm_cell.values{1,1},2)
                [v1{i}]=setdiff(perm_cell.values{cases+ii_pat,1}{i},perm_cell.values{perm_vals,1}{i});  %#ok<*NASGU,*AGROW> % identify the seeds
            end
            if ~isempty(v1{2}) && isempty(v1{1}) && isempty(v1{3})&& isempty(v1{4})&& isempty(v1{5})&& isempty(v1{6}) && isempty(v1{7})&& isempty(v1{8})&& isempty(v1{9}) && isempty(v1{10}) && isempty(v1{11}) && isempty(v1{12})
                cummulative_name=1+cummulative_name;
                accum{cummulative_name,cases}=perm_cell.OutNames{perm_vals,1}; % accumulate in the same column the names with same values but different seed
            end
        end
        
        
    end
    % loading the windfields and calculate the mean for the different seeds
    for col_names=1:size(accum,2)
        for row_names =1:size(accum,1)
            LoadWF(row_names,col_names)= load([input.LidarOutput_dir accum{row_names,col_names}]);
            TI_DATA_lidar (row_names,col_names)= LoadWF(row_names,col_names).Output.statistics.U.lidar.TI_mean;
            TI_DATA_fullWF(row_names,col_names)= LoadWF(row_names,col_names).Output.statistics.U.fullWF.TI_mean;
        end
    end
    TI_Qlundar.TI_mean_lidar  = mean (TI_DATA_lidar);
    TI_Qlundar.TI_mean_WF     = mean (TI_DATA_fullWF);
    for ind_err=1:size(TI_Qlundar.TI_mean_WF,2)
        TI_Qlundar.error (1,ind_err) = 100*(abs(TI_Qlundar.TI_mean_lidar(1,ind_err)-TI_Qlundar.TI_mean_WF(1,ind_err))/TI_Qlundar.TI_mean_WF(1,ind_err)); % error [%]
    end 
    %save data
    save_data_full_path = [input.Qlundar_TI accum{row_names,col_names} '_TIout.mat'];
    save(save_data_full_path,'TI_Qlundar');
    disp(['Turbulence intensity for QlunDAR: ' accum{row_names,col_names} ' has been processed (' datestr(datetime) '):' ])
    disp('Creating virtual lidar TI output for QlunDAR finished successfully')
    
    ii_pat=size(perm_cell.OutNames,1)/size(input.PatternNames,2); % Doesn't work if the patterns have need a smarter way of doing that!!!!!
end
