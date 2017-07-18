%=========================================================================
%
%     Nonaparametric density of ftse returns 
%
%=========================================================================

clear all;
clc;

% Load equities returns data on ftse 20 Nov. 1973 to 23 July 2001 
load ftse.mat 

yt = data;
t  = length(yt);
y  = -0.10:0.001:0.10;

% Generate the Student t distribution for comparison
mu  = mean(yt);
sig = std(yt);
nu  = 5;

f_studt = tpdf( (y - mu)/sig,nu )*(1/sig);
% tmp = ( y - mu )/sig;
% f_studt = gamma((nu+1)/2)*(1/gamma(0.5))*(1/gamma(nu/2))*(1/sqrt(nu))*(1 + tmp.^2/nu).^(-(1+nu)/2)*(1/sig);

% Estimate the nonparametric density
f = zeros(length(y),1);
h = 1.06*std(yt)*t^(-1/5);

for i = 1:length(y);
    z     = (y(i) - yt)./h;    
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
plot(yt,'-k');
title('(a) Data');
ylabel('$y_t$');
xlabel('$t$');
ylim([-0.12 0.08]);
xlim([0 7000]);
set(gca,'YTick',[-0.12 -0.08 -0.04  0.0 0.04 0.08])
set(gca,'XTick',[1000 3000 5000 7000])
box off;

%--------------------------------------------------------%
% Panel (b)
subplot(1,2,2);
plot(y,f,'--k','LineWidth',0.5);
hold on;
plot(y,f_studt,'-k','LineWidth',0.5);
hold off;
title('(b) Distributions');
ylabel('$f(y_t)$');
xlabel('$y$');
ylim([0,60]);
set(gca,'YTick',[0 10 20 30 40 50 60]);

box off;

%laprint(1,'ftseret','options','factory');


