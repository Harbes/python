%=========================================================================
%
%   Maximum likelihood estimation of the stationary distribution of the 
%   Vasciek model of interest rates using Ait Sahalia's (1996) data.
%
%=========================================================================
clear all;
clc;

% Load data (5505x4 array called eurodata, 1 Jun 1973 - 25 Feb 1995)
load eurodata.mat
rt = eurodata(:,4)*100;

% Maximum likelihood estimates of the stationary distribution 
mu_r   = mean( rt );
sig2_r = mean( (rt - mu_r).^2 );

disp( ['MLE of mean of stationary distribution:      ', num2str(mu_r)]);
disp( ['MLE of variance of stationary distribution: ', num2str(sig2_r)]);


% Compute stationary density from -5% to 35%
r     = -5:0.1:25;
fnorm = normpdf( ( r - mu_r )/sqrt(sig2_r) )/sqrt(sig2_r);
prob  = normcdf( ( 0 - mu_r )/sqrt(sig2_r) );

disp( ['Probability of a negative interest rate:     ', num2str(prob)]);

%***     Generate graph    ***%

% Switch off TeX interpreter and clear figure
% First figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

plot(r,fnorm,'-k')
ylabel('f(r)')
xlabel('Interest Rate');
set(gca,'ytick',[]);
axis tight;

% Print the tex file to the relevant directory
%laprint(1,'statdensity','options','factory');
