%=========================================================================
%
%   Program to demonstrate properties of the LM test
%
%=========================================================================
clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',12345) );

% Sample sizes
t = 200;

% Simulate the data (discard first 100 simulated observations)
mu   = 0.0;            
phi1 = 0.8;           
si1  = 0.3;
vt   = randn(t+101,1);
yt   = recserar( mu + trimr(vt,1,0) + si1*trimr(vt,0,1) , 0.0 , phi1);  
yt   = trimr(yt,100,0);   

% Demonstrate the equivalence result  
% LM test of AR(1)
% First regression
y = yt;   
x = ones(length(y),1);
v = y - x*(x\y);

% Second regression
y   = trimr(v,1,0);    
x   = [ones(length(y),1)  trimr(v,0,1) ];
e   = y - x*(x\y);
y   = bsxfun(@minus,y,mean(y));
rsq = 1 - (e'*e)/(y'*y);     
lm  = (t-1)*rsq;

disp(' ')
disp(['LM statistic (phi1 = 0.0)      = ', num2str(lm) ]);
disp(['p-value                        = ', num2str(1-chi2cdf(lm,1)) ]);
disp(' ')

% LM test of MA(1)
% First regression
y = yt; 
x = ones(length(y),1);
e = y - x*(x\y);

% Second regression
y   = trimr(e,1,0);    
x   = [ones(length(y),1)   trimr(e,0,1) ];
e   = y - x*(x\y);
y   = bsxfun(@minus,y,mean(y));
rsq = 1 - (e'*e)/(y'*y);     
lm  = (t-1)*rsq;

disp(' ')
disp(['LM statistic (si1 = 0.0)       = ', num2str(lm) ]);
disp(['p-value                        = ', num2str(1-chi2cdf(lm,1)) ]);
disp(' ')


% Demonstrate singularity result
% LM test (joint AR and MA)   
% First regression
y = yt;   
x = ones(length(y),1);
v = y - x*(x\y);

% Second regression
y = trimr(v,1,0);                   
x = [ones(length(y),1)   trimr(yt,0,1)   trimr(v,0,1) ];

disp('Correlation of explanatory variables at second stage')
disp( corrcoef( x(:,[2 3]) ) );
