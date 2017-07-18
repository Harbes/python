
%***************************************************************************
%***
%***     Program to demonstrate the biasedness property of MLE for the 
%***     normal distribution and an adjustment to achieve unbiasedness.
%***
%***    ** Note that the results will not match the numbers reported in the
%***     text exactly because of differences in the Gauss and Matlab
%***     random number generation.**
%***
%***************************************************************************

clear all;
clc;

state = 123457;
rand('state', state);
randn('state', state);


mue  = 1;       % Population mean        
sig2 = 2;       % Population variance    
r    = 20000;   % Number of replications                            
t    = 5;       % Sample size            


%     Generate sample means assuming population mean is unknown      

u   = randn(t,r);                            % Generate N(0,1) random numbers          
y   = mue + sqrt(sig2)*u;                    % Generate N(mue,sig2) random numbers     
tmp = repmat(mean(y),size(y,1),1);

vary          = sum( (y - tmp).^2 )/t;       % Compute the MLEs of the r samples       
vary_unbiased = sum( (y - tmp).^2 )/(t-1);   % Compute the unbiased estimates of the r samples       

fprintf('Results based on unknown population mean\n');
fprintf('****************************************\n');
fprintf('Theoretical value                               = %f\n', sig2);
fprintf('Simulated expectation of the MLE                = %f\n', mean(vary));
fprintf('Simulated expectation of the unbiased estimator = %f\n\n\n', mean(vary_unbiased));


%     Generate sample means assuming population mean is known       

u = randn(t,r);                         	% Generate N(0,1) random numbers         
y = mue + sqrt(sig2)*u;                     % Generate N(mue,sig2) random numbers  

vary          = sum( (y - mue).^2 )/t;     	% Compute the MLEs of the r samples      
vary_unbiased = sum( (y - mue).^2 )/(t-1); 	% Compute the unbiased estimates of the r samples      


fprintf('Results based on known population mean\n');
fprintf('**************************************\n');
fprintf('Theoretical value                               = %f\n', sig2);
fprintf('Simulated expectation of the MLE                = %f\n', mean(vary));
fprintf('Simulated expectation of the unbiased estimator = %f\n\n\n', mean(vary_unbiased));


