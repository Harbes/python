
%=======================================================================
%
%    Program to estimate a Normal model with unknown mean 
%   and known variance equal to one
%
%=======================================================================

clear all;
clc;

% Read in the data for y
y = [ 1; 2; 5; 1; 2; ];
t = length( y );   	% Define the sample size

% Compute the MLE
theta = mean(y);

disp('MLE');
disp('---');
disp(theta);


% Define the log of the likelihood at each observation

lnlt = -0.5*log(2*pi) - 0.5*(y - theta).^2;	


% Evaluate log like at each obs

disp('Log like at MLE for each obs');
disp('----------------------------');
disp([y lnlt])

% Evaluate log likelihood 
 
disp('Log Likelihood at MLE');
disp('---------------------');
disp(mean(lnlt));

 
% Evaluate the gradient at the MLE  for each obs
 
g_t = y - theta;
 
disp('Gradient of the log like at MLE for each obs');
disp('--------------------------------------------');
disp(g_t);
 
% Evaluate the gradient at the MLE
g = mean(g_t);
 
disp('Gradient of the log likelihood at MLE');
disp('-------------------------------------');
disp(g);


% Evaluate the Hessian at the MLE  for each obs
h = -t;

disp('Hessian of the log likelihood at MLE');
disp('------------------------------------');
disp(h);

