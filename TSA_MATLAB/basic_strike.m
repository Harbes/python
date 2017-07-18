%=========================================================================
%
%   Program to model the duration of strikes per annum in the U.S.
%   Data are from Kennan (1985)
%
%=========================================================================

clear all;
clc;

% Read the data: US strike data from 1968 to 1976
load strike.mat;

yt = data(:,1);                

% Maximum likelihood estimate of strike duration 
theta = mean(yt);

disp( ['MLE of mean strike duration (in days) = ' num2str(theta)]);


% Plot the estimated distribution of strike durations  
y = 0:10:300;
f = (1/theta)*exp(-y/theta);

figure(1)
plot(y,f);

% Plot the histogram 
figure(2)
hist(yt,21);

