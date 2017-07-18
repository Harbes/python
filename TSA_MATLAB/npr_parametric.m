%=========================================================================
%
%     Parametric solutions of the nonlinear example 
%
%=========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) );

% Simulate the model
t  = 500;
ut = 0.1*randn(t,1);                % N(0,0.1^2) 
xt = -2 + 4*rand(t,1);              % U(-2,2) 
mx = 0.3*exp( -4*(xt + 1).^2 ) + 0.7*exp( -16*(xt - 1).^2 );  
yt  = mx + ut;

% Estimate a linear model
y     = yt;
x     = [ones(t,1) xt];
b1    = x\y;
yfit1 = x*b1;

% Estimate a nonlinear (polynomial) model
y     = yt;
x     = [ones(t,1) xt xt.^2 x.^3 xt.^4];
b2    = x\y;
yfit2 = x*b2;

%**********************************************************************
%***
%***     Generate graphs
%***
%**********************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

tmp = sortrows( [xt mx yfit1 yfit2],1 );

%--------------------------------------------------------%
% Panel (a)
subplot(1,2,1);
plot(tmp(:,1),tmp(:,2),'-k');
hold on;
plot(tmp(:,1),tmp(:,3),'-.k','LineWidth',1.0);
hold off;
title('(a) Linear');
ylabel('$m(x_t)$');
xlabel('$x_t$');
%legend('True','Estimated','Location','NorthWest')
%legend boxoff;
box off;

%--------------------------------------------------------%
% Panel (b)
subplot(1,2,2);
plot(tmp(:,1),tmp(:,2),'-k');
hold on;
plot(tmp(:,1),tmp(:,4),'-.k','LineWidth',1.0);
title('(b) Nonlinear');
ylabel('$m(x_t)$');
xlabel('$x_t$');
%legend('True','Estimated','Location','NorthWest')
%legend boxoff
box off;

%laprint(1,'figpara','options','factory');
