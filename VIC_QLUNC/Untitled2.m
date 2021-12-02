clear all;
close all;
clc
input=[11.1442194798231,11.4749720978737,11.6699736413956,12.2677283728838,12.2984786162853,12.5249804091454];
% input=[11,11,12,11,11];

x=linspace(-length(input),length(input),size(input,2));
pdf=fitdist(input','Normal');
x_gaus=normpdf(x,0,1);

mean_weight = sum(x_gaus.*input)/sum(x_gaus)

figure,
plot (x,x_gaus)
title("inputs")














% figure,
% plot (input,gauss)
% title("gaus")


B = imgaussfilt(input,pdf.sigma)

w=gausswin(length(input))

fiat=conv(input,w)

mean_weight = sum(w.*fiat)/sum(w)

