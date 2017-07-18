%=========================================================================
%
%   Program to simulate the Poisson autoregressive model 
%   using the binomial thinning operator.
%
%=========================================================================
clear all
clc
    
RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) );

% Parameters
t   = 10;                                                    
rho = 0.3;                                 
lam = 3.5;

% Generate poisson random numbers
u = poissrnd(lam,t,1);                          

% Initialize y by choosing the median of draws from the unconditional distribution 
% to ensure that the starting value is positive)  
y = ones(t,1)*median(poissrnd(lam/(1-rho),t,1));

for i = 2:t

    % Generate y_t-1 uniform random numbers   
    e = rand(y(i-1),1);             
   
    % Sum bernoulli random draws  
    b_thin = sum( e < rho );      

    % Generate realisation of y at time t   
    y(i) = b_thin + u(i);           

    disp(['Iteration     = ',num2str(i) ]); 
    disp(['y(i-1)        = ',num2str(y(i-1)) ]);
    disp(['b_thin        = ',num2str(b_thin) ]); 
    disp(['u(i)          = ',num2str(u(i)) ]);
    disp(['y(i)          = ',num2str(y(i)) ]);
    disp( ' ' );
    

end