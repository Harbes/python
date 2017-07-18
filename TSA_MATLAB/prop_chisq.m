
%***************************************************************************
%***
%***     Monte Carlo program to demonstrate the Central Limit Theorem
%***     where 10000 numbers are drawn from a Chi-squared distribution
%***     with one degree of freedom.
%***
%***     For each sample of size 5, the sample mean is computed and the
%***     standardized random variable constructed where the population
%***     mean and variance are equal to 1 and 2 respectively for the
%***     Chi-squared distribution with one degree of freedom
%***
%***************************************************************************

clear all;
clc;
clf;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',123457) )

        
t = 5;                               	% Sample size                   
r = 1000;                           	% Number of draws                 
rnorm = randn(t,r);                 	% Generate N(0,1) random numbers   
rchi1 = rnorm.^2;                   	% Chi-squared (1) random numbers   
z = sqrt(t)*(mean(rchi1)' - 1)/sqrt(2); % Standardized random variable

hist(z,21);                         	% Plot the histogram with 21 bars


%***************************************************************************
%***                                                                        
%***     Note that the last four lines can be written on one line as:       
%***                                                                        
%***             hist( (mean(rnorm.^2)' -1 )/sqrt(2/n), 21);                
%***                                                                        
%***************************************************************************

fprintf('Replications       = %f\n', r);
fprintf('Sample size        = %f\n', t);
fprintf('Mean               = %f\n', mean(z)');
fprintf('Variance           = %f\n', std(z)'^2);
fprintf('Standard deviation = %f\n\n', std(z)');

zsort = sort(z,1);                 % Sort the data 

fprintf('The empirical distribution\n');
fprintf('**************************\n');
fprintf(' 1.0 per cent      = %f\n', zsort(0.01*r,1));
fprintf(' 2.5 per cent      = %f\n', zsort(0.025*r,1));
fprintf(' 5.0 per cent      = %f\n', zsort(0.05*r,1));
fprintf('10.0 per cent      = %f\n', zsort(0.10*r,1));
fprintf('50.0 per cent      = %f\n', zsort(0.50*r,1));
fprintf('90.0 per cent      = %f\n', zsort(0.90*r,1));
fprintf('95.0 per cent      = %f\n', zsort(0.95*r,1));
fprintf('97.5 per cent      = %f\n', zsort(0.975*r,1));
fprintf('99.0 per cent      = %f\n\n', zsort(0.99*r,1));
