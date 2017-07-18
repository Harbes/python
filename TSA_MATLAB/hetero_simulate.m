% ========================================================================
%
%    Simulating a model of heteroskedasticity 
%
% ========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',123456) )

% Simulate the model
t     = 5000;
alpha = 1;
beta  = 2;
delta = 1;
gam   = 0.5;

x = randn(t,1);                      	                    
w = (0.0:0.1:0.1*(t-1))';               %     wt is a time trend   
u = sqrt(delta + gam*w).*randn(t,1); 	                          
y = alpha + beta*x + u;             	                                     
           

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
subplot(2,2,1);
plot(x,y,'-k');
title('(a)');
xlabel('$x_t$');
ylabel('$y_t$');
box off;

%--------------------------------------------------------%
% Panel (b)
subplot(2,2,2);
plot(w,y,'-k');
title('(b)');
xlabel('$w_t$');
ylabel('$y_t$');
xlim([0,500]);
box off;

%--------------------------------------------------------%
% Panel (c)
subplot(2,2,3);
plot(w,y.^2,'-k');
title('(c)');
xlabel('$w_t$');
ylabel('$y^{2}_{t}$');
xlim([0,500]);
box off;

%--------------------------------------------------------%
% Panel (d)

bols = [ones(t,1) x]\y;
e = y - [ones(t,1) x]*bols; 

subplot(2,2,4);
plot(w,e.^2,'-k');
title('(d)');
xlabel('$x_t$');
ylabel('$e^{2}_{t}$');
xlim([0,500]);
box off;

%laprint(1,'heterosim','options','factory');

