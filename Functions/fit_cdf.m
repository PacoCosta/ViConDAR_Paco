
function Vfinal=fit_cdf(X,Y)

% MEthod 1
% fcn = @(b,x) normcdf(x, b(1), b(2));                    % Objective Function
% NRCF = @(b) norm(Y - fcn(b,X));                     % Norm Residual Cost Function
% B = fminsearch(NRCF, [0; 10]);                          % Estimate Parameters
% X = linspace(min(X), max(X));
% Y = fcn(B,X);
% [val_unique_Fvals,ind_unique_Fvals]=unique(Y);
% try
%     Vfinal = interp1(val_unique_Fvals,X(ind_unique_Fvals),0.5,'linearinterp');
% catch
%     Vfinal = 0;
% figure(1)
% plot(X, Y, 'pg')
% hold on
% plot(Xplot, fcn(B,Xplot))
% hold off
% grid
% text(-50, 0.65, sprintf('\\mu = %.1f\n\\sigma = %.1f', B))

% Method 2    
[val_unique_Fvals,ind_unique_Fvals]=unique(Y); % Take unique values to calculate the interpolation. The cdf is normalised using the maximum
fit_cdf = fit (val_unique_Fvals',X(ind_unique_Fvals)','linearinterp'); %Fit the cdf to get the value in cdf=0.5
Vfinal  = fit_cdf(.5); % Value in cdf=.5
    
end