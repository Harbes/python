%=========================================================================
%
%   Maximum likelihood estimation of the transitional distribution of the 
%   Vasciek model of interest rates using Ait Sahalia's (1996) data.
%
%=========================================================================
clear all;
clc;

% Load data (5505x4 array called eurodata, 1 Jun 1973 - 25 Feb 1995)
load eurodata.mat
rt = eurodata(:,4)*100;

% Regress r(t) on r(t-1)
y  = trimr(rt,1,0);
x  = [ones( length( y ),1 )  trimr(rt,0,1) ];

theta = x\y;
e     = y - x*theta;
sig2  = mean(e.^2);

disp( ['MLE of alpha:      ', num2str( theta(1) )]);
disp( ['MLE of rho:        ', num2str( theta(2) )]);
disp( ['MLE of beta:       ', num2str( theta(2)-1 )]);
disp( ['MLE of sigma^2:    ', num2str( sig2 )]);


disp( ' ' );

disp( ['Minimum interest rate:      ', num2str( min( rt ) )]);
disp( ['Median interest rate:       ', num2str( median( rt ) )]);
disp( ['Maximum interest rate:      ', num2str( max( rt ) )]);

disp( ' ' );

disp( ['MLE of mean based on the transitional distribution:   ', num2str( -theta(1)/(theta(2)-1) )]);
disp( ['MLE of variance based on the transitional distribution:   ', num2str( -sig2/((theta(2)-1)*(2+theta(2)-1)) )]);

% Compute transitional densities
r = 0:0.1:30;

mu_min = theta(1) + theta(2)*min(rt);
fmin   = normpdf( ( r - mu_min )/sqrt(sig2) )/sqrt(sig2);
 
mu_med = theta(1) + theta(2)*median(rt);
fmed   = normpdf( ( r - mu_med )/sqrt(sig2) )/sqrt(sig2);

mu_max = theta(1) + theta(2)*max(rt);
fmax   = normpdf( ( r - mu_max )/sqrt(sig2) )/sqrt(sig2);

%***     Generate graph ***%

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
f1=figure(1);
clf;

plot(r,fmin,'--k',r,fmed,'-k',r,fmax,':k','LineWidth',0.75);
ylabel('f(r)')
xlabel('r');
set(gca,'ytick',[ ]);
axis tight;

%legend('Minimum','Median','Maximum','Location','Best')
%legend2latex( f1 )

% Print the tex file to the relevant directory
%laprint(1,'transdensity','options','factory');

