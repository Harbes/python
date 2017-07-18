%=========================================================================
%
%   Program to compute gradient, Hessian, information matrix and OPG
%=========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',12) )
t    = 5000;         

% Normal distribution  
mu   = 1;                  
sig2 = 1;                  
u    = randn(t,1);              
y    = mu  + sqrt(sig2 )*u;     
m    = mean(y);              % maximum likelihood estimate      
g    = (y - m)/sig2;         % gradient at t    
h   = -1/sig2;               % hessian at t    
j   = g'*g;                  % opg  

disp(' Normal distribution results ')
disp(['Gradient vector    = ',num2str(mean(g)) ]);
disp(['Hessian            = ',num2str(mean(h)) ]);
disp(['Information matrix = ',num2str(-mean(h)) ]);
disp(['OPG matrix         = ',num2str(j/t) ]);
disp(' ')

% Exponential distribution 
theta = 2;               

y  = theta*gamrnd(1,1,[t 1]);   
th = 1/mean(y);             % maximum likelihood estimate
g  = 1/th - y;              % gradient at t 
h  = -1/th^2;               % hessian at t
j = g'*g;                   % opg

disp(' Exponential distribution results ')
disp(['Gradient vector    = ',num2str(mean(g)) ]);
disp(['Hessian            = ',num2str(mean(h)) ]);
disp(['Information matrix = ',num2str(-mean(h)) ]);
disp(['OPG matrix         = ',num2str(j/t) ]);
disp(' ')
