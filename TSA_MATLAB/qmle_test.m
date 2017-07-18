%========================================================================
%
%   Program to look at log(X^2) and normal.
%
%========================================================================

clear all
clc


grid = (-15:0.01:5)';

% Kernel estimate of density of log(x^2)

x = sortrows( randn(10000,1) );
f1 = ksdensity( log(x.^2),grid ); 

f2 =  (1/sqrt(2*pi)).*exp( (grid-exp(grid))/2);


% Normal approximation
ny = normpdf(grid,-1.27,sqrt(4.93));
 
plot(grid, [ny f1  f2] );
legend('N(0,1)','kernel','theory');