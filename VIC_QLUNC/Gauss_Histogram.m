gaussian = @(x, m,s) exp(-0.5*((x-m)/s).^2)/(s*sqrt(2*pi));
mu = 5;                             % Arbitrary mean
sigma = 1;                          % arbitrary standard deviation
x = mu + sigma*randn(1,400);   % Arbitrary random normal data
h = histogram(x);
dh = h.BinWidth;
lo = min(x);
hi = max(x);
dxx = (hi-lo)/100;
xx = linspace(lo,hi,101);
pdf = gaussian(xx, mu, sigma);
scalefactor = sum(h.Values*dh)/(trapz(pdf)*dxx);
pdf = scalefactor*pdf;
hold on
plot(xx,pdf)

xx(1)-xx(2)
xx(34)-xx(35)
xx(78)-xx(79)


figure(1),
hold on,
plot(distan,gaussian_w,'-b'),
% histfit(gaussian_w)
grid off
hold off


x = [randn(100,1); 4+randn(50,1)];
[hts,ctrs] = hist(x)
bar(ctrs,hts,'hist')
area = sum(hts) * (ctrs(2)-ctrs(1))
xx = linspace(-3,7);
hold on; plot(xx,area*normpdf(xx,mean(x),std(x)),'r-')
f = ksdensity(x,xx);
plot(xx,area*f,'g-')
hold off

