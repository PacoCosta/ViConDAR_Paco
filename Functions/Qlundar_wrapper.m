%% Header
%
% Qlundar wrapper. Gets the names from perm_cell and calculates the mean TI
% and mean RMSE
% from the different seeds. 
%
% F.Costa
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2021



function TI_Qlundar = Qlundar_wrapper(input,perm_cell)
ii_pat=0;
for ind_pattern = 1: length(input.PatternNames)
    
    %     Obtain the names
    n_cases=(length(input.freeInp{1,2})*length(input.freeInp{3,2})*length(input.freeInp{4,2})*length(input.timestep_pat_vec {1})*length(input.timeStep_Measurements{1})) ;
    for cases = 1:n_cases
        cummulative_name = 1;
        accum{1,cases} = perm_cell.OutNames{cases+ii_pat,1};
        for perm_vals = 1:size(perm_cell.values,1)
            for i = 1:size (perm_cell.values{1,1},2)
                [v1{i}] = setdiff(perm_cell.values{cases+ii_pat,1}{i},perm_cell.values{perm_vals,1}{i});  %#ok<*NASGU,*AGROW> % identify the seeds
            end
            if ~isempty(v1{2}) && isempty(v1{1}) && isempty(v1{3})&& isempty(v1{4})&& isempty(v1{5})&& isempty(v1{6}) && isempty(v1{7})&& isempty(v1{8})&& isempty(v1{9}) && isempty(v1{10}) && isempty(v1{11}) && isempty(v1{12})
                cummulative_name = 1+cummulative_name;
                accum{cummulative_name,cases} = perm_cell.OutNames{perm_vals,1}; % accumulate in the same column the names with same values but different seed
            end
        end        
    end
    % loading the windfields and calculate the mean for the different seeds
    for col_names = 1:size(accum,2)
        for row_names = 1:size(accum,1)
            LoadWF(row_names,col_names) = load([input.LidarOutput_dir accum{row_names,col_names}]);
            % u-component
            TI_DATA_lidar_U (row_names,col_names) = LoadWF(row_names,col_names).Output.statistics.U.lidar.TI_mean;
            TI_DATA_fullWF_U(row_names,col_names) = LoadWF(row_names,col_names).Output.statistics.U.fullWF.TI_mean;
            RMSE_U (row_names,col_names)          = mean(LoadWF (row_names,col_names).Output.statistics.U.lidar.RMSE); % we mean between the different points in the pattern for the same seed wind field
            
            % v-component
            TI_DATA_lidar_V (row_names,col_names) = LoadWF(row_names,col_names).Output.statistics.V.lidar.TI_mean;
            TI_DATA_fullWF_V(row_names,col_names) = LoadWF(row_names,col_names).Output.statistics.V.fullWF.TI_mean;
            RMSE_V (row_names,col_names)          = mean(LoadWF (row_names,col_names).Output.statistics.V.lidar.RMSE);
            
            % w-component
            TI_DATA_lidar_W (row_names,col_names) = LoadWF(row_names,col_names).Output.statistics.W.lidar.TI_mean;
            TI_DATA_fullWF_W(row_names,col_names) = LoadWF(row_names,col_names).Output.statistics.W.fullWF.TI_mean;
            RMSE_V (row_names,col_names)          = mean(LoadWF(row_names,col_names).Output.statistics.W.lidar.RMSE);
            
        end
    end
    % mean of windfiled characteristics from different seeds of the same wind field
    TI_mean_lidar_U = mean (TI_DATA_lidar_U);
    TI_mean_WF_U    = mean (TI_DATA_fullWF_U);
    RMSE_mean_U     = mean (RMSE_U);
    for ind_err = 1:size(TI_mean_WF_U,2)
        error_U (1,ind_err) = 100*(abs(TI_mean_lidar_U(1,ind_err)-TI_mean_WF_U(1,ind_err))/TI_mean_WF_U(1,ind_err)); % error [%]
    end 
    %save data
    for in_save=1:size(accum,2)
        save_data_full_path        = [input.Qlundar_TI accum{1,in_save} '_TIout.mat'];
        TI_Qlundar.TI_mean_lidar_U = TI_mean_lidar_U (1,in_save);
        TI_Qlundar.TI_mean_WF_U    = TI_mean_WF_U(1,in_save);
        TI_Qlundar.error_U         = error_U(1,in_save);      % relative error in TI estimations 
        TI_Qlundar.RMSE_mean_U     = RMSE_mean_U(1,in_save);  
        
        save(save_data_full_path,'TI_Qlundar');
        disp(['Turbulence intensity for QlunDAR: ' accum{1,in_save} ' has been processed (' datestr(datetime) '):' ])
        
    end
    disp(['Creating virtual lidar TI output for QlunDAR finished successfully'])
    ii_pat=size(perm_cell.OutNames,1)/size(input.PatternNames,2); % Doesn't work if the patterns have need a smarter way of doing that!!!!!
end
