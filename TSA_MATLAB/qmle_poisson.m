%========================================================================
%
%   Program to compute the QMLE standard errors where the true 
%   distribution is Poisson and the misspecified distribution is normal
%
%========================================================================
clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',12) )


mu = 26;         
t  = 100;         

%  Generate data from the true model (Poisson)  
y = poissrnd(mu,t,1);                  


% Estimate the true model (Bernoulli) 
theta0 = mean(y);

% Gradient of the true model (at each observation)  
g0 = y/theta0 - 1;

% Hessian of the true model                             
h0 = mean( -y/theta0^2 );

% Information of the true model                       
i0 = 1/theta0;


% OPG of the true model 
j0 = g0'*g0/t;

% Estimate the misspecified model (normal)    
theta1  = mean(y);

% Gradient of the misspecified model (at each observation)     
g1 = y - theta1;

% Hessian of the misspecified model                       
h1 = -1;

% Information of the misspecified model   
i1 = 1;


% OPG of the misspecified model (independence)    
j1 = g1'*g1/t;


disp(['Covariance matrix true (i)            = ' num2str(inv(i0)/t) ]);
disp(['Covariance matrix true (h)            = ' num2str(inv(-h0)/t) ]);
disp(['Covariance matrix true (j)            = ' num2str(inv(j0)/t) ]);
disp(' ');
disp(['Covariance matrix misspecified (i)    = ' num2str(inv(i1)/t) ]);
disp(['Covariance matrix misspecified (h)    = ' num2str(inv(-h1)/t) ]);
disp(['Covariance matrix misspecified (j)    = ' num2str(inv(j1)/t) ]);
disp(['Covariance matrix misspecified (qmle) = ' num2str(inv(i1)*j1*inv(i1)/t) ]);
