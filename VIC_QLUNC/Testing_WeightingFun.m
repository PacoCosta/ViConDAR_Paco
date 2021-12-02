% close all 
% clear all
% clc
% fwhm_zr   = 2*30;
fwhm  = 2*input.distance_av_space;   
% focus   = 150;
% offset = 100;
% interval_of_confidence=1e4;
interval_of_confidence = 4*sigmas;
% distan = linspace(focus-offset,focus+offset,1e4);
% distan = linspace(-50,50,accuracy);
% sigma=fwhm_zr/2.355;
sigmas = fwhm/(2*sqrt(2*log(2)));
distan     = linspace(-interval_of_confidence,interval_of_confidence,size(VFinalTotal_TimeInt3,1));
gaussian = (1/(sigmas*sqrt(2*pi)))*exp(-0.5*((distan-focus)/sigmas).^2);

% Find the half max value.
halfMax = (min(gaussian) + max(gaussian)) / 2;
% Find where the data first drops below half the max.
index1 = find(gaussian >= halfMax, 1, 'first');

% Find where the data last rises above half the max.
index2 = find(gaussian >= halfMax, 1, 'last');
fwhm_ind = index2-index1; % FWHM in indexes.
fwhmx = distan(index2) - distan(index1)

%cumulative probability:
cumulative_probability = cumsum(gaussian)*(distan(3)-distan(2));
sum_probability        = sum(gaussian)*(distan(2)-distan(1))

%% analysis:
% velocity vector:
a = 10;
b = 12;
v_vec = (b-a).*rand(1,accuracy) + a;
% v_vec= sin(0.1:.01:100);
% v_vec=10*ones(1,accuracy);
V_mean=sum(v_vec)/length(v_vec)
% Gaussian factors, velocities and distances above the fwhm
[~,ind] = find(gaussian >= halfMax);
probe_distance=distan(ind);

v_vec2=v_vec(ind);
gaussian2=gaussian(ind);

% take some values for the analysis
%%%%%%%%%%%%%%%%%%%%%%%%
% filter=randi([0 1], 1,length(probe_distance2));
% filter =ones(length(probe_distance2),1)';
% probe_distance_filt=nonzeros(probe_distance2.*filter);
% v_vec_filt=nonzeros(v_vec.*filter);
% gaussian_filt=nonzeros(gaussian2.*filter);
%%%%%%%%%%%%%%%%%%%%%%%%%%

for step=1:length(gaussian2)
    filt_vec{1,step}(:,:)=ones(1,length(probe_distance));
    filter_step{1,step}(:,:)=1:step:length(probe_distance);
    filt_vec{1,step}(filter_step{1,step})=0;
    
    v_vec_filt{1,step}=filt_vec{1,step}.*v_vec2;
    gaussian_filt{1,step}=filt_vec{1,step}.*gaussian2;
    VFinalTotal(:,step) = sum(gaussian_filt{1,step}.*v_vec_filt{1,step})/sum(gaussian_filt{1,step});
end


%% plotting

% Plot analysis
figure(1),hold on,
plot(1:length(VFinalTotal),VFinalTotal)
plot(1:length(VFinalTotal),repelem(mean(VFinalTotal,'omitnan'),length(VFinalTotal)),'r-','Linewidth',1.5)
plot(1:length(VFinalTotal),repelem(mean(VFinalTotal+0.00005*VFinalTotal,'omitnan'),length(VFinalTotal)),'r--',1:length(VFinalTotal),repelem(mean(VFinalTotal-0.00005*VFinalTotal,'omitnan'),length(VFinalTotal)),'r--','Linewidth',2)
plot(1:length(VFinalTotal),repelem(mean(VFinalTotal+0.0001*VFinalTotal,'omitnan'),length(VFinalTotal)),'g--',1:length(VFinalTotal),repelem(mean(VFinalTotal-0.0001*VFinalTotal,'omitnan'),length(VFinalTotal)),'g--','Linewidth',2)

title ('Velocity value convergence Vs number of points in the probe volume')
xlabel ('n° points within the probe volume [-]')
ylabel ('U vel [m/s]')
set(gca,'FontSize',23)
hold off

% plot gaussian
figure(2),
hold on,
plot(distan,gaussian,'-b'),
plot(distan,linspace(halfMax,halfMax,length(distan)))
grid off
hold off

% Plot cummulative probability
figure(3), plot(distan,cumulative_probability);




