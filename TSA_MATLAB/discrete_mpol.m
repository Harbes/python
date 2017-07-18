%=========================================================================
%
%   Plot US Federal Funds target rate
%
%=========================================================================

clear all
clc

% Read the data: loads a time series object "target"
% US weekly yields starting 1st week of Feb 1984 -- first week of June 1997

load FedFunds

dates = getabstime( target );


%**********************************************************************
%***
%***     Generate graph
%***
%**********************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

plot(datenum(dates),target.Data,'-k','LineWidth',0.75);
datetick('x','yyyy')
ylabel('Federal Funds Rate');
xlabel('Time')
axis( [-Inf, Inf,2,12] );
box off


% Print the tex file to the relevant directory
%laprint(1,'fedrate','options','factory');

  
