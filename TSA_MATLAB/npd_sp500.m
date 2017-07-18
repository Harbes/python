%=========================================================================
%
%     Estimated distribution of S&P 500
%
%=========================================================================

clear all;
clc;

% Load equities data on 2&p500 9 Feb. 1954 to 23 July 2001    
load sp500.mat; 

yt = sp500;
t  = length(yt);

% Generate the Student t distribution for comparison    
mu      = mean(yt);
sig     = std(yt);
nu      = 5;
y       = (-0.25:0.001:0.25)';
z       = ( y - mu )/sig;
f_studt = gamma((nu+1)/2)*(1/gamma(0.5))*(1/gamma(nu/2))*(1/sqrt(nu))*(1 + z.^2/nu).^(-(1+nu)/2)*(1/sig);

% Estimate the nonparametric density     
f = zeros(length(y),1);
h  = 1.06*std(yt)*t^(-1/5);
for i = 1:length(y);
    z     = ((y(i) - yt)./h);    
    f(i) = mean( normpdf(z)./h );
end

%**********************************************************************
%***
%***     Generate graph
%***
%**********************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

%--------------------------------------------------------%
% Panel (a)
subplot(1,2,1);
plot(yt,'-.k','LineWidth',1.0);
ylabel('$y]]t[[$');
xlabel('$t$');
axis tight
box off;

%--------------------------------------------------------%
% Panel (b)
subplot(1,2,2);
plot(y,f,'-.k','LineWidth',1.0);
hold on;
plot(y,f_studt,'-k','LineWidth',0.5);
hold off;
ylabel('$f(y)$');
xlabel('$y$');
axis tight
box off;
