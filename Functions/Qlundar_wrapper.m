%% Header
%
% Qlundar wrapper. Gets the names fron perm_cell and calculates the TI
%
% F.Costa
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2021



function [TI_Qlundar]=Qlundar_wrapper(input,perm_cell)

% Obtain the names
% n_cases=(length(input.freeInp{1,2})*length(input.freeInp{3,2})*length(input.freeInp{4,2})*length(input.timestep_pat_vec {1})*length(input.timeStep_Measurements{1})) ;
% for first=1:n_cases
%     cummulative_name=1;
%     accum{1,first}=perm_cell.OutNames{first,1};
%     for second=1:size(perm_cell.values,1)
%         for i = 1:size (perm_cell.values{1,1},2)
%             [a{i},v{i}]=(setdiff(perm_cell.values{first,1}{i},perm_cell.values{second,1}{i}));  %#ok<*NASGU,*AGROW> % identify the seeds
%         end
%         if ~isempty(a{2}) && isempty(a{1}) && isempty(a{3})&& isempty(a{4})&& isempty(a{5})&& isempty(a{6}) && isempty(a{7})&& isempty(a{8})&& isempty(a{9}) && isempty(a{10}) && isempty(a{11}) && isempty(a{12})
%             cummulative_name=1+cummulative_name;
%             accum{cummulative_name,first}=perm_cell.OutNames{second,1}; % accumulate in the same column the names with same values but different seed
%         end
%     end
% end
%
%
% % loading the windfields and calculate the mean for the different seeds
% for col_names=1:size(accum,2)
%     for row_names =1:size(accum,1)
%         LoadNam{row_names,col_names}= load([input.LidarOutput_dir accum{row_names,col_names}]);
%         TI_DATA_lidar {row_names,col_names}= LoadNam{row_names,col_names}.Output.statistics.U.lidar.TI_mean;
%         TI_DATA_fullWF{row_names,col_names}= LoadNam{row_names,col_names}.Output.statistics.U.fullWF.TI_mean;
%     end
%     TI_Qlundar.TI_mean_lidar  = mean (TI_DATA_lidar {row_names,col_names});
%     TI_Qlundar.TI_mean_WF     = mean (TI_DATA_fullWF{row_names,col_names});
%     TI_Qlundar.error          = (1-(TI_Qlundar.TI_mean_lidar/TI_Qlundar.TI_mean_WF));
%     %save
%
%     save_data_full_path = [input.Qlundar_TI accum{row_names,col_names} '_TIout.mat'];
%     save(save_data_full_path,'TI_Qlundar');
%     disp(['Turbulence intensity for QlunDAR: ' accum{row_names,col_names} ' has been processed (' datestr(datetime) '):' ])
%     disp('Creating virtual lidar TI output for QlunDAR finished successfully')
% end
indTI=1;
for i_TI=input.freeInp{4,2}
    
    i_TI=num2str(i_TI);
    if size(i_TI,2)<2
        i_TI=['0' i_TI];
    end
    indTp=1;
    for i_Tp =input.timestep_pat_vec{1}
        
        i_Tp=num2str(i_Tp);
        if size(i_Tp,2)<2
            i_Tp=['0' i_Tp];
        end
        
        i_Tp = strrep(i_Tp,'.','d');
        
        ind=1;
        for i_Dav = {'06' ,'30','60','90','120','150d0','180','210','240','270','300'} % user puts here the number of points appearing in the name ('Dav')
            
            load([input.Qlundar_TI 'QlunDAR_Sh00_SD10_V15_TI' i_TI '_1P_Single_Tp' i_Tp '_Tm00_Fd250_DAv' i_Dav{1} '_TIout.mat']);
            TI_lidar2(:,ind)=(TI_Qlundar.TI_mean_lidar);
            TI_fullWF2(:,ind)=(TI_Qlundar.TI_mean_WF);
            ind=ind+1;
        end
        
        TI_lidar1{indTp}=TI_lidar2;
        TI_fullWF1{indTp}=TI_fullWF2;
        
        indTp=indTp+1;
    end
    TI_lidar{indTI}=TI_lidar1;
    TI_fullWF{indTI}=TI_fullWF1;
    indTI=1+indTI;
end
sj=0;
