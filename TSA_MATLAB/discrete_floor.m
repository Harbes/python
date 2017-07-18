%=========================================================================
%
%   Money demand equation with a floor interest rate
%
%=========================================================================
clear all
clc

RandStream.setGlobalStream( RandStream('mt19937ar','seed',42) );


t  = 20;
x  = seqa(0,1,t)';
b0 = 10;
b1 = -0.5;

m = b0 + b1*x;
s = 1;
c = 4;

y  = m + s*randn(t,1);
f1 = normcdf( (c - m)/s );
f2 = normpdf( (c - m)/s );
pred = m + s*f2./(1 - f1);


%**********************************************************************
%***
%***     Generate graph
%***
%**********************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

plot(x,m,'--k',x,c*ones(t,1),':k',x,pred,'-k')
ylabel('Interest rate (%)');
xlabel('Money');
box off

% Print the tex file to the relevant directory
%laprint(1,'floor','options','factory');



