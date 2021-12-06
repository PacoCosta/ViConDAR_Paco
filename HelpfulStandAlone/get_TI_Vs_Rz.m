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
Zr = [120   150   180   210   240   270   300];
Pattern_name        = '1P_Single';
Sufix               = '_TIout';
V                   = '15';
FD                  = '250';
Sh                  ='143';
SD                  = '01';
TI                  = [2,6,10,14,18,22,26,30];
Tp                  = [0.1,1,10];
distance_av         = [10,12.5,15,17.5,20,22.5,25];
ColorSet = varycolor(8);
ind2=1;
for i_Tp = Tp
    ind_leg=1;
    i_Tp = num2str(i_Tp);
    if size(i_Tp,2)<2
        i_Tp = ['0' i_Tp];
    end
    i_Tp = strrep(i_Tp,'.','d');
    ind0=1;
    indColor=1;
    for i_TI=TI
        i_TI = num2str(i_TI); %#ok<*FXSET>
        if size(i_TI,2)<2
            i_TI = ['0' i_TI]; %#ok<*AGROW>
        end
        ind1=1;
        for i_Zr = Zr % user puts here the number of points appearing in the name ('Dav')
            i_str_Zr=num2str(i_Zr);
            load([Qlunda_TI_directory nameBase '_Sh' Sh '_SD' SD '_V' V '_TI' i_TI '_' Pattern_name '_Tp' i_Tp '_Tm00_Fd' FD '_DAv' i_str_Zr Sufix '.mat']);
            Output_TI_Qlundar.TI_lidar_U{ind0}(:,ind1)   = TI_Qlundar.TI_mean_lidar_U;
            Output_TI_Qlundar.TI_fullWF_U{ind0}(:,ind1)  = TI_Qlundar.TI_mean_WF_U;
            Output_TI_Qlundar.TI_error_U{ind0}(:,ind1)   = TI_Qlundar.error_U;
            Output_TI_Qlundar.RMSE_mean_U {ind0}(:,ind1) = TI_Qlundar.RMSE_mean_U;
            ind1=ind1+1;
        end
        ind0=ind0+1;
        legCell{1,ind_leg}=['e(TI) - TI = ' i_TI , ' Tp = ' num2str(Tp(ind2))];
        legCell2{1,ind_leg}=['RMSE(V) - TI = ' i_TI , ' Tp = ' num2str(Tp(ind2))];
        ind_leg=ind_leg+1;
        save([save_dir nameBase '_TI' i_TI '_Tp' i_Tp ], 'Output_TI_Qlundar')
    end
%     ind2=ind2+1;
    % Plotting
    figure,hold on
    for i = [1:size( Output_TI_Qlundar.TI_lidar_U,2)]
        plot(distance_av,Output_TI_Qlundar.TI_error_U{i},'-','color', ColorSet(indColor,:))
        
        title ('e(TI) and RMSE(U) in lidar estimates')
        grid on
        set(gca,'Fontsize',19)
        legend (legCell,'Interpreter','None','FontSize',10)
        xlabel('Rayleigh distance [m]')
        
        yyaxis right
        ylabel('[-]')
        
        plot(distance_av,Output_TI_Qlundar.RMSE_mean_U{i},'--','color', ColorSet(indColor,:))
        legend (legCell2,'Interpreter','None','FontSize',10)
        
        yyaxis left
        ylabel('[%]')
        indColor=indColor+1;
    end
    hold off
    %contours
    % Performing matrices for the contours
    for ind_cont=1:8
        matCont_RMSE{ind2}(ind_cont,:)=Output_TI_Qlundar.RMSE_mean_U{ind_cont};
        mat_Cont_errTI{ind2}(ind_cont,:)=Output_TI_Qlundar.TI_error_U{ind_cont};
    end
    
    % RMSE contour
    CLimits_RMSE = [floor(min(min(matCont_RMSE{ind2}))), round(max(max(matCont_RMSE{ind2})),4)];
    figure, contourf(distance_av,TI,matCont_RMSE{ind2});
    title(['RMSE (U) - Tp = ' num2str(Tp(ind2) )])
    c0=colorbar;
    c0.Label.String = 'RMSE [-]';
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
