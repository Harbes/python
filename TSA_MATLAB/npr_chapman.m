%=========================================================================
%
%    Program to reproduce the results in Chapman and Pearson, 
%    Journal of Finance (2000), p.355-388. 
%
%=========================================================================
clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',12345) );

t = 7500;

kappa = 0.21459;
theta = 0.085711;
sigma = 0.07830;
delta = 1/250;

% Simulate the interest rate 
r = theta + zeros(t,1);

for i = 2:t

   r(i) = r(i-1) + kappa*(theta - r(i-1))*delta + sigma*sqrt(r(i-1))*sqrt(delta)*randn(1,1);

end

% Nonparametric regression of the mean
dr = r(2:end) - r(1:end-1);
yt = dr/delta;
xt = r(1:end-1);
h  = 0.023;
x  = (0.0:0.002:0.2)';

fx = zeros(length(x),1);
fxy = zeros(length(x),1);
for i = 1:length(x);
    z      = ((x(i) - xt)./h);    
    fx(i)  = mean( normpdf(z)./h );
    fxy(i) = mean( normpdf(z).*yt./h );
end
mx  = fxy./fx;

yt_mean   = yt;
mx_mean   = mx;
true_mean = kappa*(theta - x);

% Nonparametric regression of the variance 
dr = r(2:end) - r(1:end-1);
yt = dr.^2/delta;
xt = r(1:end-1);
h  = 0.023;

fx = zeros(length(x),1);
fxy = zeros(length(x),1);
for i = 1:length(x);
    z      = ((x(i) - xt)./h);    
    fx(i)  = mean( normpdf(z)./h );
    fxy(i) = mean( normpdf(z).*yt./h );
end
mx  = fxy./fx;

yt_var   = yt;
mx_var   = mx;
true_var = sigma^2*x;

%**************************************************************************
%**
%**     Generate graphs
%**
%**************************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

%--------------------------------------------------------%
% Panel (a)
subplot(2,2,1);
plot(yt_mean,'-k');
title('(a) $y_t = (r_t - r_{t-1})/\Delta$');
xlabel('$t$');
ylabel('$y_{t}$');
axis tight
box 'off'

%--------------------------------------------------------%
% Panel (b)
subplot(2,2,2);
plot(x,true_mean,'-k', ...
     x,mx_mean,'-.k');
title('(b) Estimated Mean');
ylabel('$m(x_{t})$');
xlabel('$x_{t}$');
axis tight
box 'off'

%--------------------------------------------------------%
% Panel (c)
subplot(2,2,3);
plot(yt_var,'-k');
title('(c) $y_t = (r_t - r_{t-1})^2/\Delta$');
xlabel('$t$');
ylabel('$y_t$');
axis tight
box 'off'

%--------------------------------------------------------%
% Panel (d)
subplot(2,2,4);
plot(x,true_var,'-k', ...
     x,mx_var,'-.k');
title('(c) Estimated Variance');
xlabel('$x_{t}$');
ylabel('$m(x_{t})$');
axis tight
box 'off'

% laprint(1,'chapman','options','factory');
