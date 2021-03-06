%% Header
% Qlundar wrapper. Obtain TI Vs Zr vector to plot TI differences according to
% the probe length
% F.Costa
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2021

% function get_TI_Vs_Rz(Zr)
close all
clear all
clc
%read from:
Qlunda_TI_directory = '..\Qlundar_TI\';
%Save in:
save_dir            = '..\Output_Qlundar_TI\';
nameBase            = 'Qlundar';
% Introduce wf inputs (base on the name of the wf stored in '..\Qlundar_TI\')
Zr = [120];%[120,150,180,210,240,270,300]; % user puts here the number distance appearing in the name ('Dav'=Zr*TruncationDistance)
Pattern_name        = 'RHI';
Sufix               = '_TIout';
V                   = '15';
FD                  = '150';
Sh                  ='143';
SD                  = '01';
TI                  = [2,30];%[2,6,10,14,18,22,26,30];
Tp                  = [1];%[0.1,1,10];
distance_av         = [10];%[10,12.5,15,17.5,20,22.5,25]; 
ColorSet            = varycolor(size(TI,2));
ind2=1;
for i_Tp = Tp
    ind_leg=1;
    i_Tp = num2str(i_Tp);
    if size(i_Tp,2)<2
        i_Tp = ['0' i_Tp];
    end
    i_Tp = strrep(i_Tp,'.','d');
    ind0=1;
    
    for i_TI=TI
        i_TI = num2str(i_TI); %#ok<*FXSET>
        if size(i_TI,2)<2
            i_TI = ['0' i_TI]; %#ok<*AGROW>
        end
        ind1=1;
        for i_Zr = Zr 
            
            i_str_Zr=num2str(i_Zr);
            if contains(i_str_Zr,'.')
                i_str_Zr = strrep(i_str_Zr,'.','d');
            end
            load([Qlunda_TI_directory nameBase '_Sh' Sh '_SD' SD '_V' V '_TI' i_TI '_' Pattern_name '_Tp' i_Tp '_Tm00_Fd' FD '_DAv' i_str_Zr Sufix '.mat']);
            Output_TI_Qlundar.TI_lidar_U{ind0}(:,ind1)   = TI_Qlundar.TI_mean_lidar_U;
            Output_TI_Qlundar.TI_fullWF_U{ind0}(:,ind1)  = TI_Qlundar.TI_mean_WF_U;
            Output_TI_Qlundar.TI_error_U{ind0}(:,ind1)   = TI_Qlundar.error_U;
            Output_TI_Qlundar.RMSE_mean_U {ind0}(:,ind1) = TI_Qlundar.RMSE_mean_U;
            ind1=ind1+1;
        end
        ind0=ind0+1;
        legCell{1,ind_leg}=[[char(949),'(TI)'] ' - TI = ' i_TI ];
        legCell2{1,ind_leg}=['RMSE(u) - TI = ' i_TI ];
        ind_leg=ind_leg+1;
        save([save_dir nameBase '_TI' i_TI '_Tp' i_Tp ], 'Output_TI_Qlundar')
    end
    %     ind2=ind2+1;
    % Plotting
    figure,hold on,grid on
    indColor=1;
    for i = 1:size( Output_TI_Qlundar.TI_lidar_U,2)
        plot(distance_av,Output_TI_Qlundar.TI_error_U{i},'o','color', ColorSet(indColor,:));
        title ([char(949),'(TI)' , ' (Tp = ' num2str(Tp(ind2)),'s)']);
        
        set(gca,'Fontsize',19)
        xlabel('Rayleigh distance [m]')
        ylabel('[%]')
        indColor=indColor+1;
    end
    legend(legCell,'Fontsize',10)
    
    % Plotting RMSE velocity
    figure, hold on,grid on   
    indColor=1;
    for i = 1:size( Output_TI_Qlundar.TI_lidar_U,2)         
        plot(distance_av,Output_TI_Qlundar.RMSE_mean_U{i},'o','color', ColorSet(indColor,:));       
        set(gca,'Fontsize',19)
        title (['RMSE(u)' , ' (Tp = ' num2str(Tp(ind2)),'s)']);
        ylabel('[m/s]')
        xlabel('Rayleigh distance [m]')
        indColor=indColor+1;
        legend()
    end
    legend(legCell2,'Fontsize',10)
   
    %contours
    % Performing matrices for the contours
    for ind_cont=1:size(TI,1)
        matCont_RMSE{ind2}(ind_cont,:)=Output_TI_Qlundar.RMSE_mean_U{ind_cont};
        mat_Cont_errTI{ind2}(ind_cont,:)=Output_TI_Qlundar.TI_error_U{ind_cont};
    end
    
    % RMSE contour
%     CLimits_RMSE = [floor(min(min(matCont_RMSE{ind2}))), round(max(max(matCont_RMSE{ind2})),4)];
    CLimits_RMSE = [0, 1.65];
    figure, contourf(distance_av,TI,matCont_RMSE{ind2});
    title(['RMSE (U) - Tp = ' num2str(Tp(ind2) )])
    c0=colorbar;
    c0.Label.String = 'RMSE [m/s]';
    xlabel('Rayleigh distance [m]')
    ylabel('TI[%]')
    set(gca,'BoxStyle','full','CLim',CLimits_RMSE,'Layer','top','Fontsize',19);
    
    % error in TI contour
    CLimits_eTI = [floor(min(min(mat_Cont_errTI{ind2}))), round(max(max(mat_Cont_errTI{ind2})),4)];
    figure, contourf(distance_av,TI,mat_Cont_errTI{ind2});
    title(['e(TI) [%] - Tp = ' num2str(Tp(ind2) )])
    c1=colorbar;
    c1.Label.String = 'e[%]';
    xlabel('Rayleigh distance [m]')
    ylabel('TI[%]')
    set(gca,'BoxStyle','full','CLim',CLimits_eTI,'Layer','top','Fontsize',19);
    ind2=ind2+1;
end