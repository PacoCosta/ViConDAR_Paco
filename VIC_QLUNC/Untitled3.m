close all 
clear all
clc
fwhm_zr   = 2*14.23;
focus   = 250;
offset = 50;
distan = linspace(-offset+focus-fwhm_zr,offset+focus+fwhm_zr,1e4);

sigma=fwhm_zr/2.355;
gaussian = (1/(sigma*sqrt(2*pi)))*exp(-0.5*((distan-focus)/sigma).^2);

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
v_vec = (b-a).*rand(1,1e4) + a;
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

figure,hold on,plot(1:length(VFinalTotal),VFinalTotal)
plot(1:length(VFinalTotal),repelem(mean(VFinalTotal,'omitnan'),length(VFinalTotal)))
hold off
%% plot gaussian

figure,
hold on,
plot(distan,gaussian,'-b'),
plot(distan,linspace(halfMax,halfMax,length(distan)))
grid on
hold off
figure, plot(distan,cumulative_probability);