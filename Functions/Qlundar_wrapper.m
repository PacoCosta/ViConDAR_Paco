% function TI_err=Qlundar_wrapper(input,perm_cell,)

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
% end




% end