% ========================================================================
%
%    Simulate Brownian motion and compare continuous and discrete data
%
% ========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',12345) );

dt = 1/60;                      % Minute data
n  = 240;                       % 10 days = 240hrs 
t = n/dt;                       % Sample size

% Parameters
mu   = 0.0;                     % Mean
sig2 = 1.0;                     % Variance

% Wiener process
dw = sqrt(dt)*randn(t,1);                                        

% Generate continous data (minutes)
y(1) = 0.0;

for i = 1:t-1
    
    y(i+1)=y(i)+mu*dt + sqrt(sig2)*dw(i);
end

y10 = y(10:10:t);               % Generate 10 minute data by choosing every 10th observation
yhr = y(10:60:t);               % Generate hourly data by choosing every 60th observation 
ydy = y(1440:1440:t);           % Generate daily data by choosing every 1440th observation         **/


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
ind = (1:1:t)/1440;
subplot(2,2,1);
plot(ind,y,'-k');
title('(a) Minute Data');
ylabel('$y_t$');
xlabel('$t$ days');
set(gca,'XTick',0:1:10);
set(gca,'YTick',-10:10:30);
ylim([-10,30]);
xlim([0,10]);
box off;


%--------------------------------------------------------%
% Panel (b)
subplot(2,2,2);
ind = (10:10:t)/1440;
plot(ind,y10,'-k')
title('(b) Ten Minute Data');
ylabel('$y_t$');
xlabel('$t$ days');
set(gca,'XTick',0:1:10);
set(gca,'YTick',-10:10:30);
ylim([-10,30]);
xlim([0,10]);
box off;

%--------------------------------------------------------%
% Panel (c)
ind = (10:60:t)/1440;
subplot(2,2,3);
plot(ind,yhr,'-k')
title('(c) Hourly Data');
ylabel('$y_t$');
xlabel('$t$ days');
set(gca,'XTick',0:1:10);
set(gca,'YTick',-10:10:30);
ylim([-10,30]);
xlim([0,10]);
box off;

%--------------------------------------------------------%
% Panel (d)
ind = (1440:1440:t)/1440;
subplot(2,2,4);
plot(ind,ydy,'-k')
title('(c) Daily Data');
ylabel('$y_t$');
xlabel('$t$ days');
set(gca,'XTick',0:1:10);
set(gca,'YTick',-10:10:30);
ylim([-10,30]);
xlim([0,10]);
box off;

%laprint(1,'brown','options','factory');
