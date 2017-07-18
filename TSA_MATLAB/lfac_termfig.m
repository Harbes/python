%=========================================================================
%
%   Program to plot term structure data 
%
%=========================================================================

clear all
clc

% Read the data:
% US daily zero coupon yields starting 10/04/1988 and ending 12/28/2001
% 	The variables are:
%
%		1.	tcm3m
%		2.	tcm1y
%		3.	tcm3y
%		4.	tcm5y
%		5.	tcm7y
%		6.	tcm10y


load lfac_usdata.mat

data = usdata;


% Set up dates
t = 1988+10/12:1/250:2002.057333333333;
t = t(1:length(data));



%**********************************************************************
%***
%***     Generate graphs
%***
%**********************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

%--------------------------------------------------------%
% Panel (a)
subplot(3,2,1)
plot(t,data(:,1),'-k','LineWidth',1);
title('(a) 3 Month Yield');
ylabel('\%');
xlabel('t');
xlim( [ 2 10 ] )
box off
axis tight


%--------------------------------------------------------%
% Panel (b)
subplot(3,2,2)
plot(t,data(:,2),'-k','LineWidth',1);
title('(b) 1 Year Yield');
ylabel('\%');
xlabel('t');
box off
axis tight


%--------------------------------------------------------%
% Panel (c)
subplot(3,2,3)
plot(t,data(:,3),'-k','LineWidth',1);
title('(c) 3 Year Yield');
ylabel('\%');
xlabel('t');
box off
axis tight


%--------------------------------------------------------%
% Panel (d)
subplot(3,2,4)
plot(t,data(:,4),'-k','LineWidth',1);
title('(d) 5 Year Yield');
ylabel('\%');
xlabel('t');
box off
axis tight


%--------------------------------------------------------%
% Panel (e)
subplot(3,2,5)
plot(t,data(:,5),'-k');
title('(e) 7 Year Yield','LineWidth',1);
ylabel('\%');
xlabel('t');
box off
axis tight


%--------------------------------------------------------%
% Panel (f)
subplot(3,2,6)
plot(t,data(:,6),'-k');
title('(f) 10 Year Yield','LineWidth',1);
ylabel('\%');
xlabel('t');
xlim([ 2 10] );
%xTick( [2 4 6 8]);
box off
axis tight

        
        
