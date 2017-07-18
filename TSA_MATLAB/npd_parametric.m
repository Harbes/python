%=========================================================================
%
%     Normal distribution example using different parametric distributions
%
%=========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',1) );

% Simulate the data from a normal distribution       
t   = 200;        %      Sample size     
mu  = 0;
sig = 3;
yt  = mu + sig*randn(t,1);

% Generate the population (normal) distribution       
y      = (-15:0.1:15)';
f_norm = normpdf( ( y - mu )/sig )/sig;

% Estimate the parametric density assuming normality        
m       = mean(yt);
s       = std(yt);
f_npara = normpdf( ( y - m )/s )/s;

% Estimate the parametric density assuming Student t       
m       = mean(yt);
s       = std(yt);
v       = 2*s^2/(s^2-1);
f_spara = (1+((y-m)/s).^2/v).^(-(v+1)/2)*gamma((v+1)/2)/(sqrt(pi*v)*gamma(v/2)*s);

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
plot(yt,'-.k','LineWidth',1.0);
title('(a) Data');
ylabel('$y(t)$');
xlabel('$t$');
box off;

%--------------------------------------------------------%
% Panel (b)
subplot(1,2,2);
plot(y,f_norm,'-k','LineWidth',0.5);

title('(b) Distribution');
ylabel('$f(y)$');
xlabel('$y$');
box off;


% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(2);
clf;

%--------------------------------------------------------%
% Panel (a)
subplot(1,2,1);
plot(y,f_npara,'-.k','LineWidth',1.0);
hold on;
plot(y,f_norm,'-k','LineWidth',0.5);
hold off;
title('(a)');
ylabel('$f(y)$');
xlabel('$y$');
box off;

%--------------------------------------------------------%
% Panel (b)
subplot(1,2,2);
plot(y,f_spara,'-.k','LineWidth',1.0);
hold on;
plot(y,f_norm','-k','LineWidth',0.5);
hold off;
title('(b)');
ylabel('$y(t)$');
xlabel('$t$');
box off;