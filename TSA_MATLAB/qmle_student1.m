%=========================================================================
%
%   Program to demonstrate the effects of missspecifying the 
%   log-likelihood function on the relationship between the 
%   information matrix and the outer product of gradients.
%
%   True model is student t, misspecified model is normal.
%
%=========================================================================

    
clear all
clc
   
RandStream.setDefaultStream( RandStream('mt19937ar','seed',12) )

    
mu  = 1;     % Population mean      
sig = 1;     % Population standard deviation  
gam = 5;    % Degrees of freedom
t   = 500000;

% Generate data from the true model (Student t)
v = tinv( rand(t,1),gam );           
y = mu + sig*sqrt( (gam-2)/gam )*v;

% Estimate the misspecified model (normal)    
m  = mean(y);
s2 = mean( (y - m).^2 );

% Compute gradients of the misspecified model
g1 = (y - m)/s2;
g2 = -0.5/s2 + 0.5*(y - m).^2/s2^2;
g = [ g1  g2 ] ;

 
% Compute Hessian and information matrix of the misspecified model  
H = zeros(2,2);

H(1,1) = -1/s2;
H(1,2) = -mean(y - m)/s2^2;
H(2,1) = H(1,2);
H(2,2) = 0.5/s2^2 - mean((y - m).^2)/s2^3;
    
I  = -H;
    
% Compute opg of the misspecified model (normal) 
J = g'*g/t;
    
disp('Information Matrix');
disp('------------------');
disp( I  );
    
  
disp('OPG Matrix' )
disp('------------------');
disp( J );
    

disp('Difference' )
disp('------------------');
disp( (I-J) );

 
