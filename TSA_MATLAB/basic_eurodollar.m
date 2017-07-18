%=========================================================================
%
%   Interest rate application - daily Eurodollar rates (1973-1995)
%   Data used by Ait-Sahalia (RFS 1996)
%
%=========================================================================

clear all;
clc;

% Load data
load eurodata.mat

%***********************************************************************
%***
%***     Generate graph
%***
%***********************************************************************

% Switch off TeX interpreter and clear figure
% First figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

plot(eurodata(:,3),eurodata(:,4)*100,'-k')
ylabel('$\%$');
xlabel('$t$');
axis tight;
set(gca,'ytick',[ 4 8 12 16 20 24]);

% Print the tex file to the relevant directory
laprint(1,'eurodollar','options','factory');

