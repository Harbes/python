
%***************************************************************************
%***
%***     Program to find the MLEs using the Newton-Raphson algorithm
%***
%***************************************************************************

clear all;
clc;

state = 123;
randn('state', state);


% Part (a): Simulate the model 

x = [ 1; 2; 4; 5; 8; ];

beta = 1.0;
sig2 = 4.0;
t = size(x,1);
y = beta*x + sqrt(sig2)*randn(t,1);

% (3) - Evaluate gradient and Hessian at theta0 

beta  = 1.0;
sig2  = 4.0;
theta = beta | sig2;

g       = zeros(2,1);
g(1,1)  = sum( (y - beta*x).*x )/sig2;
g(2,1)  = -0.5*t/sig2 + 0.5*sum( (y - beta*x).^2 )/sig2^2;

h       = zeros(2,2);
h(1,1)  = -sum(x.^2)/sig2;
h(1,2)  = -sum( (y - beta*x).*x )/sig2^2;
h(2,1)  = -sum( (y - beta*x).*x )/sig2^2;
h(2,2)  = 0.5*t/sig2^2 - sum( (y - beta*x).^2 )/sig2^3;

% Average log-likelihood at theta0

a_0     = -0.5*log(2*pi) - 0.5*log(sig2) - 0.5*mean( (y - beta*x).*2 )/sig2;

% Newton-Raphson

while g'*g > 0.001

theta = theta - inv(h)*g;
beta  = theta(1);
sig2  = theta(2);

g       = zeros(2,1);
g(1,1)  = sum( (y - beta*x).*x )/sig2;
g(2,1)  = -0.5*t/sig2 + 0.5*sum( (y - beta*x).^2 )/sig2^2;

h       = zeros(2,2);
h(1,1)  = -sum(x.^2)/sig2;
h(1,2)  = -sum( (y - beta*x).*x )/sig2^2;
h(2,1)  = -sum( (y - beta*x).*x )/sig2^2;
h(2,2)  = 0.5*t/sig2^2 - sum( (y - beta*x).^2 )/sig2^3;

end

% Optimum Average log-likelihood at final theta

a_1     = -0.5*log(2*pi) - 0.5*log(sig2) - 0.5*mean( (y - beta*x).*2 )/sig2;

fprintf(1,'Average log-likelihood at theta0 = %f\n\n', a_0);
fprintf(1,'After applying newton raphson\n');
fprintf(1,'Final value of beta              = %f\n', beta);
fprintf(1,'Final value of variance          = %f\n', sig2);
fprintf(1,'Final average log-likelihood     = %f\n', a_1);




