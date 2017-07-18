%=========================================================================
%
%   Program to demonstrate the properties of an artificial neural network
%
%=========================================================================

clear all;
clc;

RandStream.setGlobalStream( RandStream('mt19937ar','seed',123457) );

t = 200;      

% Set the parameters of the ANN      
phi    = 0.4;           
gam    = 2.0;
delta0 = -2.0;
delta1 = 2.0;


% Create predictions  						     
ylag =-4:0.1:4;

f = 1./( 1+exp(-( delta0 + delta1*ylag)) );

ylin  = phi*ylag;
ynlin = gam*f;

yhat = ylin + ynlin;



%**********************************************************************
%***
%***     Generate graph
%***
%**********************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

plot(ylag,yhat,'-k',  ...
     ylag,ylin,'--k', ...
     ylag,ynlin,':k');

ylabel('$y_t$');
xlabel('$y_{t-1}$');
set(gca,'YTick',[-2 -0 2 4]);
set(gca,'XTick',[-4 -2 0 2 4]);
box off


% Print the tex file to the relevant directory
%laprint(1,'annfig','options','factory');
