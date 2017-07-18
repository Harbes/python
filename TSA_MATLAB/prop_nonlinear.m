%=========================================================================
%
%   Program to demonstrate multiple roots of a regression model
%   with nonlinear parameterisation
%
%=========================================================================
clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',12) )

t = 100;                            

% Parameters
beta0 = 2.5;
beta1 = -1.5;
sig2  = 1.0;

% Generate data
u = randn(t,1);                       
x = rand(t,2);                    
y = beta0 + beta1*x(:,1) + beta1^2*x(:,2) + sqrt(sig2)*u;


ngrid = 299;
lnl   = zeros(ngrid,1);                 
g     = zeros(ngrid,2);                 
b1    = seqa(-2.1,0.01,ngrid);         


for i = 1:ngrid
    
    b0 = beta0;
    z      = (y - b0 - b1(i)*x(:,1) - b1(i)^2*x(:,2))/sig2;
    lnl(i) = - log(2*pi) - 0.5*log(sig2) - 0.5*mean(z.^2);
    g(i,1) = mean( (y - b0 - b1(i)*x(:,1) - b1(i)^2*x(:,2)) );
    g(i,2) = mean( (y - b0 - b1(i)*x(:,1) - b1(i)^2*x(:,2)).*(x(:,1) + 2*b1(i)*x(:,2)) );

end



figure(1)
subplot(1,2,1)
plot(b1,[g(:,2) zeros(ngrid,1)]);
title('(a) Gradient');
ylabel('$G_T(\theta)$')
xlabel('$\theta$')

subplot(1,2,2)
plot(b1,lnl)
title('(b) Log-likelihood');
ylabel('$L_T(\theta)$');
xlabel('$\theta$')
