% ========================================================================
%
%     Cauchy example
%
% ========================================================================


% Simulate the data from a Cauchy distribution       


clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) );

t = 100;         %      Sample size            
r = 1000;        %      Number of replications 

nue = 1;
yt  = 1./(0.5 - atan(rand(t,r))/pi);  % Inverse of cauchy cdf. 

tstat = sqrt(t)*mean(yt)./std(yt);

% Generate the asymptotic (normal) distribution       
y      = -5:0.01:4.99;
f_norm = normpdf(y,0,1);

% Estimate the nonparametric density       
h  = 1.06*std(tstat')*t^(-1/5);
wt = normpdf(((y - tstat)./h))'./h;
f  = mean(wt); 

%**************************************************************************
%**
%**     Generate graph
%**
%**************************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

%--------------------------------------------------------%
% Panel (a)
subplot(1,2,1);
plot(tstat,'-.k','LineWidth',1.0);

title('(a) T-stats');
ylabel('lt');
xlabel('lr');
%axis([-5 5 0 0.5]);
%set(gca,'XTick',[-4 -2 0 2 4 ]);
%set(gca,'YTick',[0.0 0.1 0.2 0.3 0.4 0.5]);
box off;

%--------------------------------------------------------%
% Panel (b)
subplot(1,2,2);
plot(y,f,'-.k','LineWidth',1.0);
hold on;
plot(y,f_norm,'-r','LineWidth',1.0);
title('(b) Distributions');
ylabel('f(t)');
xlabel('t');
%axis([-5 5 0 0.5]);
%set(gca,'XTick',[-4 -2 0 2 4 ]);
%set(gca,'YTick',[0.0 0.1 0.2 0.3 0.4 0.5]);
box off;



