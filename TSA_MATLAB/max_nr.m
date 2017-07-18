
%***************************************************************************
%***
%***     Program to find the MLEs using one iteration of Newton-Raphson
%***
%***************************************************************************

clear all;
clc;

state = 123;
rand('twister', state);
randn('state', state);


% (1) - Simulate the model 

x = [ 1; 2; 4; 5; 8 ];

beta = 1.0;
sig2 = 4.0;
t = size(x,1);
y = beta*x + sqrt(sig2)*randn(t,1);


% (3) - Evaluate gradient and Hessian at theta0 

beta_0  = 1.0;
sig2_0  = 4.0;
theta_0 = beta_0 | sig2_0;

g       = zeros(2,1);
g(1,1)  = sum( (y - beta_0*x).*x )/sig2_0;
g(2,1)  = -0.5*t/sig2_0 + 0.5*sum( (y - beta_0*x).^2 )/sig2_0^2;

h       = zeros(2,2);
h(1,1)  = -sum(x.^2)/sig2_0;
h(1,2)  = -sum( (y - beta_0*x).*x )/sig2_0^2;
h(2,1)  = -sum( (y - beta_0*x).*x )/sig2_0^2;
h(2,2)  = 0.5*t/sig2_0^2 - sum( (y - beta_0*x).^2 )/sig2_0^3;


% (4) - Average log-likelihood at theta0

a_0     = -0.5*log(2*pi) - 0.5*log(sig2_0) - 0.5*mean( (y - beta_0*x).*2 )/sig2_0;

% (5) - Newton-Raphson update

theta_1 = theta_0 - inv(h)*g;
beta_1  = theta_1(1);
sig2_1  = theta_1(2);

% (6) - Average log-likelihood at theta1

a_1     = -0.5*log(2*pi) - 0.5*log(sig2_1) - 0.5*mean( (y - beta_1*x).*2 )/sig2_1;


fprintf(1,'Average log-likelihood at theta0 = %f\n', a_0);
fprintf(1,'Average log-likelihood at theta1 = %f\n\n', a_1);


