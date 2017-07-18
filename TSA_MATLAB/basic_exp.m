%=========================================================================
%
%   Program to estimate an exponential model and plot the 
%   log-likelihood function.
%
%=========================================================================

clear all;
clc;

% Read in the data for y  
y = [ 2.1, 2.2, 3.1, 1.6, 2.5, 0.5];

% Alternatively read in data
% load exp.mat

t     = length(y);           
theta = 1/mean(y);
lnlt  = log(theta) - theta.*y;
gt    = 1/theta - y;
ht    = -1/theta^2;

disp(['Sum of y  = ', num2str(sum(y))]);
disp(['Mean of y = ', num2str(mean(y))]);
disp(['MLE       = ', num2str(theta)]);
disp(['Log-like  = ', num2str(mean(lnlt))]);
disp(['Gradient  = ', num2str(mean(gt))]);
disp(['Hessian   = ', num2str(mean(ht))]);

disp('    yt        lnlt      gt        ht');
disp('   -------------------------------------');
disp([y' lnlt' gt' ones(t,1)*ht]);




% *************************************************************************
% ***
% ***     Generate graph
% ***
% *************************************************************************

theta = 0.001:0.01:2.5;
lnl   = t*log(theta) - theta.*sum(y);

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;


subplot(1,2,1)
plot(theta,lnl,'-k');
title('(a) Log-likelihood function');
set(gca,'YTick',-40:5:-10);
%set(gca,'XTick',0:0.5:2.5);
% ytics(-40,0.0,10,0);
% xtics(0,3,0.5,0);

ylabel('$\ln L (\theta)$');
xlabel('$\theta$');
%axis tight;
set(gca,'LineWidth',1)
box off;

subplot(1,2,2)
plot(theta,exp(lnl)*100000,'-k');

% Labels
title('(b) Likelihood function');
ylabel('$L (\theta) \times 10^5$');
xlabel('$\theta$');

% Axes
set(gca,'LineWidth',1)
%axis tight;
box off;
ylim([0,4]);
set(gca,'YTick',0.5:0.5:4.0);
%set(gca,'XTick',0:0.5:2.5);

%laprint(1,'explike','options','factory');
