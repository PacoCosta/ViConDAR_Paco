% Header


function VFinalTotal_Time = V_spec(input,VFinalTotal_TimeInt3,WeightFun,distan)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing method:
% binranges=[floor(min(VFinalTotal_TimeInt3)):.1:ceil(max(VFinalTotal_TimeInt3))];
% % binranges=[0:.1:20];
% mean_VFinalTotal_TimeInt5=mean(VFinalTotal_TimeInt3,'omitnan');
% max_VFinalTotal_TimeInt5=max(VFinalTotal_TimeInt3,'omitnan');
%
% [bincounts,counts]=histc(VFinalTotal_TimeInt3,binranges);
% figure,
% bar(binranges,bincounts,'histc')
%
% xline(mean_VFinalTotal_TimeInt5,'--g')
% xline(mode(VFinal),'--r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate the histogram:


% binwidth=.1;
% floorS1 = floorS(min(VFinalTotal_TimeInt3),2); % floor the lower limit with specific decimal, here 3 decimals
% ceilS1  = ceilS(max(VFinalTotal_TimeInt3),2);% ceil the upper limit with specific decimal, here 3 decimals
%
% n_bins=length(floorS1:binwidth:ceilS1);
% heights=histcounts(VFinalTotal_TimeInt3.*WeightFun,n_bins);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
binwidth=.1;
i_bin=1;

floorS1 = floorS(min(VFinalTotal_TimeInt3),2); % floor the lower limit with specific decimal, here 3 decimals
ceilS1  = ceilS(max(VFinalTotal_TimeInt3),2);% ceil the upper limit with specific decimal, here 3 decimals

binrang = floorS1:binwidth:ceilS1; %13:binwidth:16;%
% t=floorS1:binwidth:ceilS1;
n=length(binrang);
t = linspace(binrang(1)-binwidth/2,binrang(end)+binwidth/2,n);
dt=diff(t);
dt(end+1)=dt(1);

for ind_bin = floorS1:binwidth:ceilS1 %13:binwidth:16 %
    %     binrang(i_bin)=ind_bin;
    iv=1;
    for iii = 1:length(VFinalTotal_TimeInt3)
        if (ind_bin<VFinalTotal_TimeInt3(iii)) && (VFinalTotal_TimeInt3(iii) <=ind_bin+binwidth)
            %             DSpectrum(iv,i_bin)      = VFinalTotal_TimeInt3(iii)*WeightFun(iii)*abs(distan(iii));
            DSpectrum(iv,i_bin)      = WeightFun(iii)*abs(distan(iii));
            
        else
            %         Vel_sprectrum(iv,i_bin) = 0;
            DSpectrum(iv,i_bin)      = 0;
        end
        iv=iv+1;
    end
    i_bin=i_bin+1;
end
SumSpectrum=sum(DSpectrum);
Fvals = cumsum(SumSpectrum.*dt);
F = spline(binrang+binwidth/2,Fvals);
DF = fnder(F);  % computes its first derivative

% [counts,V_levels]=histcounts(nonzeros(VFinalTotal_TimeInt3'.*WeightFun),'BinWidth',bin_width);
% [counts1,V_levels1]=histcounts(nonzeros(Spectrum),'BinWidth',bin_width);

%% Median, centroid and max of the histogram

for ind_peak_method=1:size(input.peak_detection_method,2)
    if strcmpi(input.peak_detection_method{ind_peak_method},"median")
        %Median
        [g,b]=unique(Fvals/max(Fvals));
        
        ff=fit (g',t(b)','linearinterp');
        VFinalTotal_Time = ff(.5);
    elseif strcmpi (input.peak_detection_method{ind_peak_method},"maximum")
        %Max
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dd=fnplt(DF, 'g', 2);
        [~,n_max]=(max(dd(2,:)));
        VFinalTotal_Time=dd(1,n_max);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    elseif strcmpi (input.peak_detection_method{ind_peak_method},"centroid")
        %Centroid
        VFinalTotal_Time=sum(SumSpectrum.*binrang)/sum(SumSpectrum);
        
    elseif strcmpi(input.peak_detection_method{ind_peak_method},"mean")
        %Mean
        VFinalTotal_Time = sum(WeightFun'.*VFinalTotal_TimeInt3)/sum(WeightFun); %#ok<*AGROW>
    end
end

%% Plot the histograms
% figure,
% f=bar(binrang,Hist_size_Spec,'histc');

% xline(meanBinned_WindSpeed_Vel,'--g','LineWidth',1.5,'Displayname','mean')
% xline(maxBinned_WindSpeed_Vel,'--r','LineWidth',1.5,'Displayname','max')
% % xline(centroidBinned_WindSpeed_Spec,'--k','LineWidth',1.5,'Displayname','cen')
% title('Velocity spectra','Fontsize',25)
% xlabel('Wind speed [m/s]','Fontsize',21);
% ylabel('counts [-]','Fontsize',21);
% legend
%
% histogram(VFinalTotal_TimeInt3,'BinWidth',0.1)
% xline(meanBinned_WindSpeed_Vel,'--g','LineWidth',2,'Displayname','mean');
% xline(maxBinned_WindSpeed_Vel,'--r','LineWidth',2,'Displayname','max');
% % xline(centroidBinned_WindSpeed_Spec,'--k','LineWidth',1.5,'Displayname','cen')
% title('Velocity spectra','Fontsize',25);
% xlabel('Wind speed [m/s]','Fontsize',21);
% ylabel('Counts [-]','Fontsize',21);
% legend


% % xline(centroidBinned_WindSpeedSum_Spec,'--k','LineWidth',1.5,'Displayname','cen')
% title('Doppler spectra','Fontsize',25);
% xlabel('Wind speed [m/s]','Fontsize',21);
% ylabel('Signal strenght [-]','Fontsize',21);
% legend
% d=0;
% [~,maxVert]=max(f.Vertices);
% maxBinned_WindSpeed= (f.Vertices(maxVert(2),1)+f.Vertices(maxVert(2)+1,1))/2;
% maxBinned_WindSpeed= (binrang .* Hist_size) / (Hist_size);







