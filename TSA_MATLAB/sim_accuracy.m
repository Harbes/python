% ========================================================================
%
%   Estiamte AR(1) coefficient by simulation of a MA(1) model
%
% ========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) );


theta = 0.5;        
t     = 250;        
h     = 1;         
n     = t*h;        

% Simulate the data for a MA(1) model    
u  = randn(n,1); 
ys = trimr(u,1,0) - theta*trimr(u,0,1);

% Estimate the first order autocorrelation coefficient    
ys   = ys - mean(ys);
rhos = trimr(ys,0,1)\trimr(ys,1,0);

disp(' ')
disp(['Sample size                        = ', num2str(t) ]);
disp(['h                                  = ', num2str(h) ]);
disp(['True population parameter (theta)  = ', num2str(theta) ]);
disp(' ')
disp(['True population parameter (rho)    = ', num2str(-theta/(1+theta^2)) ]);
disp(['Simulation estimate (rho)          = ', num2str(rhos) ]);


