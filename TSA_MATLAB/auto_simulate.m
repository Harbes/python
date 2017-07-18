%==========================================================================
%
%   Simulating a regression model with autocorrelation 
%
%==========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) )

t = 200;

% Model parameters     
beta0  = 2;
beta1  = 1;
rho1   = 0.95;
delta1 = 0.95;
sigma  = 3;

% Generate the exogenous variable
x = 0.5*(1:1:t+100)' + randn(t+100,1);        	   

%     Simulate a regression model with an AR(1) disturbance term      
v = sigma*randn(t+100,1);                                                                 
u = zeros(t+100,1);                                                                           
y = zeros(t+100,1);                          	                                                

for i=2:t+100 
    u(i) = rho1*u(i-1) + v(i);                                                                    
    y(i) = beta0 + beta1*x(i) + u(i);                                                         
end
y_ar1 = y;

% Simulate a regression model with a MA(1) disturbance term      
v = sigma*randn(t+100,1);                                                                  
u = zeros(t+100,1);                                                                              
y = zeros(t+100,1);                                                                              

for i=2:t+100
    u(i) = v(i) + delta1*v(i-1);                                                                 
    y(i) = beta0 + beta1*x(i) + u(i);                                                             
end
y_ma1 = y;

% Trim data to overcome startup problems      
y_ar1 = y_ar1(101:end,:);
y_ma1 = y_ma1(101:end,:);
x     = x(101:end,:);                      
mu    = beta0 + beta1*x;                                               

%**********************************************************************
%***
%***     Plot the series
%***
%**********************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;


%--------------------------------------------------------%
% Panel (a)
subplot(1,2,1);
plot(x,mu,'-k',...
     x,y_ar1,'.k','MarkerSize',6);
title('(a) AR(1) Regression Model');
ylabel('$y_t$');
xlabel('$x_t$');
set(gca,'XTick',40:20:160);
set(gca,'YTick',20:20:160);
xlim([40,160]);
ylim([20,160]);
%set(gca,'LineWidth',1);
box off;

%--------------------------------------------------------%
% Panel (b)
subplot(1,2,2);
plot(x,mu,'-k',...
     x,y_ma1,'.k','MarkerSize',6);
title('(b) MA(1) Regression Model');
ylabel('$y_t$');
xlabel('$x_t$');
set(gca,'XTick',40:20:160);
set(gca,'YTick',20:20:160);
xlim([40,160]);
ylim([20,160]);
%set(gca,'LineWidth',1);
box off;

%laprint(1,'figauto','options','factory');


% Simulate a regression model with an AR(2) disturbance term      
t = 200;                                                                                          
beta0 = 2;
beta1 = 1;
rho1  = 0.1;
rho2  = -0.9;
sigma = 3;

x = 0.5*(1:1:t+100)' + randn(t+100,1);        	    
v = sigma*randn(t+100,1);                                                               

u0 = [0;0];
rho = [rho1;rho2];
for i = 3:length(v+1)   
    u0(i,1) = v(i-1) + rho(1)*u0(i-1) + rho(2)*u0(i-2);
end
u     = u0;
y     = beta0 + beta1*x + u;
y_ar2 = y(101:end,:);

% Simulate a regression model with an ARMA(2,2) disturbance term      
t = 200;                                                                                           
beta0  = 2;
beta1  = 1;
rho1   = 0.1;
rho2   = -0.9;
delta1 = 0.3;
delta2 = 0.2;
sigma  = 3;

x = 0.5*(1:1:t+100)' + randn(t+100,1);        	% xt is generated from a trend with normal additive errors    
v = sigma*randn(t+100,1);                       % vt is N(0,sigma^2)                                          

tmp = [0; 0; (v(3:end,:) + delta1*v(2:(end-1),:) + delta2*v(1:(end-2),:))];
u0 = [0;0];
rho = [rho1;rho2];
for i = 3:length(tmp+1)   
    u0(i,1) = tmp(i-1) + rho(1)*u0(i-1) + rho(2)*u0(i-2);
end
u      = u0;
y      = beta0 + beta1*x + u;
y_arma = y(101:end,:);

 
 
 
 