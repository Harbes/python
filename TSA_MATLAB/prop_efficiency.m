
%****************************************************************************************
%***
%***  Program to demonstrate the efficiency property of MLE for the normal distribution
%***
%***    ** Note that the results will not match the numbers reported in the
%***     text exactly because of differences in the Gauss and Matlab
%***     random number generation.**
%***
%****************************************************************************************

clear all;
clc;

state = 123457;
rand('state', state);
randn('state', state);

mue  = 1;       % Population mean         
sig2 = 2;       % Population variance     
r    = 10000;   % Number of replications                            
t    = 100;     % Sample size             


%     Generate sample means     

u = randn(t,r);             % Generate N(0,1) random numbers     
y = mue + sqrt(sig2)*u;     % Generate N(mue,sig2) random numbers    
meany   = mean(y);         	% Compute the means of the r samples     
mediany = median(y);

var_meany   = mean( (meany   - mue).^2 );
var_mediany = mean( (mediany - mue).^2 );

fprintf('Theoretical variance of the sample mean   = %f\n', sig2/t);
fprintf('Simulated variance of the sample mean     = %f\n\n', var_meany);

fprintf('Theoretical variance of the sample median = %f\n', pi*sig2/(2*t));
fprintf('Simulated variance of the sample median   = %f\n\n', var_mediany);


