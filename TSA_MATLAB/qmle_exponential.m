%========================================================================
%
%      True model is exponential, misspecified model is normal.
%
%========================================================================
    
clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234) )

t  = 500000;
mu = 0.5;                                          

% Generate data from the true model (Exponential)   
z = rand(t,1);                             
y = -mu.*log(1 - z);                       


% Estimate the misspecified model (Normal)   
m = mean(y);


% Compute gradients of the misspecified model (Normal)
g = (y - m);


disp( ['Gradients ', num2str(sum(g))] );
disp( ' ' );


% Compute hessian of the misspecified model (Normal) 
h = -t;

disp( ['Average Negative Hessian ', num2str(-h/t) ] );
disp( ' ' );

% Compute opg of the misspecified model (Normal)        
disp( ['Average OPG Matrix ', num2str(g'*g/t)] );

 

