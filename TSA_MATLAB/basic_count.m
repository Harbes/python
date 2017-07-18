%=========================================================================
%
%   Program to model the number of strikes per annum in the U.S.
%   Data are from Kennan (1985)
%
%=========================================================================

clear all;
clc;

% Number of strike per annum, 1968 to 1976  
yt = [8; 6; 11; 3; 3; 2; 18; 2;	9];       

% Maximum likelihood estimate of the number of strikes 
theta = mean(yt);
disp( ['MLE of mean number of strikes (in years) = ' num2str(theta)]);

% Plot the estimated distribution of strike numbers  
y = 0:1:20;                          
f = poisspdf(y,theta);                   

figure(1)
bar(y,f);

% Plot the histogram
figure(2) 
hist(yt,11);
