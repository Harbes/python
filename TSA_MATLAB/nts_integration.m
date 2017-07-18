%==========================================================================
%
%   	Simulate stochastic trends with different orders of integration
%
%==========================================================================
clear all
clc

RandStream.setDefaultStream( RandStream('mt19937ar','seed',15) );

% Generate simulate variables
t     = 200;
delta = 0.5;    % 0, 0.1 or 0.5 

y_noise = delta + randn(t,1);                                %      White noise
y_i0    = recserar(delta + randn(t,1),0.0,0.5);              %      Stationary process (I(0))
y_i1    = recserar(delta + randn(t,1),0.0,1.0);              %      Nonstationary process (I(1))
y_i2    = recserar(delta + randn(t,1),[0.0; 0.0],[2.0;-1.0]);%      Nonstationary process (I(2))


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
subplot(2,2,1)
plot(seqa(1,1,t),y_noise);
title('White noise');

subplot(2,2,2)
plot(seqa(1,1,t),y_i0);
title('Stationary: I(0)');


subplot(2,2,3)
plot(seqa(1,1,t),y_i1);
title('Nonstationary: I(1)');

subplot(2,2,4)
plot(seqa(1,1,t),y_i2);
title('Nonstationary: I(2)');

