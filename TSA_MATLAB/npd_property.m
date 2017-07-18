%=========================================================================
%
%     Demonstration of bias and variance of nonparametric estimators 
%     using different bandwidths. 
%
%=========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) )


% Simulate the data from a normal distribution 
t   = 500;        
mu  = 0;
sig = 1;
yt  = mu + sig*randn(t,1);

% Generate the population (normal) distribution 
y      = -5:0.1:5;
f_norm = normpdf( ( y - mu )/sig )/sig;

% Estimate the nonparametric density for alternative bandwidths
f1 = zeros(length(y),1);
h  = 1;

for i = 1:length(y);
    z     = ((y(i) - yt)./h);    
    f1(i) = mean( normpdf(z)./h );
end


f2 = zeros(length(y),1);
h = 0.1;
for i = 1:length(y);
    z     = ((y(i) - yt)./h);    
    f2(i) = mean( normpdf(z)./h );
end

plot(y,[f1 f2]);

%**********************************************************************
%***
%***     Generate graphs
%***
%**********************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

%--------------------------------------------------------%
% Panel (a)
subplot(1,2,1);
plot(y,f1,'-.k','LineWidth',1.0);
hold on;
plot(y,f_norm','-k','LineWidth',0.5);
hold off;
title('(a) Larger Bandwidth');
ylabel('$f(y_t)$');
xlabel('$y$');
axis([-5 5 0 0.5]);
set(gca,'XTick',[-4 -2 0 2 4 ]);
set(gca,'YTick',[0.0 0.1 0.2 0.3 0.4 0.5]);
box off;

%--------------------------------------------------------%
% Panel (b)
subplot(1,2,2);
plot(y,f2,'-.k','LineWidth',1.0);
hold on;
plot(y,f_norm','-k','LineWidth',0.5);
hold off;
title('(b) Smaller Bandwidth');
ylabel('$f(y_t)$');
xlabel('$y$');
axis([-5 5 0 0.5]);
set(gca,'XTick',[-4 -2 0 2 4 ]);
set(gca,'YTick',[0.0 0.1 0.2 0.3 0.4 0.5]);
box off;

%laprint(1,'m_c10figprop','options','factory');

