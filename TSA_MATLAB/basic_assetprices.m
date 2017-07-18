%=========================================================================
%
%   Program to compute the maximum likelihood estimates
%   of the asset price model.
%
%   The data consist of the Australian, Singapore and NASDAQ stock
%   market indexes for the period 3 January 1989 to 31 December 2009.
%
%=========================================================================

clear all;
clc;

% Load data
load assetprices.mat

pt = pt_aust;     % Choose the asset price 
pt = pt/1000;

% MLE based on the log-normal distribution of the share index 
% Transitions are from log price at t-1 to log price at t
y     = trimr(log(pt),1,0) - trimr(log(pt),0,1);
alpha = mean(y);
sig2  = mean( (y - alpha).^2 );

disp(['MLE of alpha based on the log-normal distribution = ', num2str(alpha)]);
disp(['MLE of sig2 based on the log-normal distribution  = ', num2str(sig2)]);

