%% Header
%
% Qlundar wrapper. Order TI for each time step of the pattern
% from the different seeds
%
% F.Costa
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2021


function get_TI_Vs_Rz(Zr)

Pattern_name = '1P_Single';
Sufix= '_TIout';
V = '15';
FD = '250';
Sh='143';
SD='01';
indTI=1;


for i_TI=input.freeInp{4,2}
    
    i_TI = num2str(i_TI);
    if size(i_TI,2)<2
        i_TI = ['0' i_TI];
    end
    indTp=1;
    for i_Tp = input.timestep_pat_vec{1}
        
        i_Tp = num2str(i_Tp);
        if size(i_Tp,2)<2
            i_Tp = ['0' i_Tp];
        end
        
        i_Tp = strrep(i_Tp,'.','d');
        
        ind=1;
        for i_Zr = Zr % user puts here the number of points appearing in the name ('Dav')
            
            load([input.Qlundar_TI input.nameBase '_Sh' Sh '_SD' SD '_V' V '_TI' i_TI '_' Pattern_name '_Tp' i_Tp '_Tm00_Fd' FD '_DAv' i_Zr{1} Sufix '.mat']);
            TI_lidar2(:,ind)  = TI_Qlundar.TI_mean_lidar;
            TI_fullWF2(:,ind) = TI_Qlundar.TI_mean_WF;
            TI_error2(:,ind)   = TI_Qlundar.error;
            ind=ind+1;
        end
        
        TI_lidar1{indTp}  = TI_lidar2;
        TI_fullWF1{indTp} = TI_fullWF2;
        TI_error1{indTp}   = TI_error2;
        indTp=indTp+1;
    end
    TI_lidar{indTI}  = TI_lidar1;
    TI_fullWF{indTI} = TI_fullWF1;
    TI_error{indTI}   = TI_error1;
    indTI = 1+indTI;
end


