%=========================================================================
%
%     Program to demonstrate the consistency property of MLE for the
%     mean of the normal distribution.
%
%=========================================================================
clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) )

t    = 500;    
mu   = 1;               % Population mean
sig2 = 2;               % Population variance  

% Generate sample means from sample of size t=1,2,3,...,t        

yb  = zeros(t,1);
for t = 1:t
    y     = mu + sqrt(sig2)*randn(t,1);
    yb(t) = mean(y);
end

%**************************************************************************
%**     
%**     Generate graph     
%**                                         
%**************************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

plot(1:1:t,mu.*ones(t,1),'-k',1:1:t,yb,'-k');
ylabel('$\bar{y}$');
xlabel('T');
box off;
axis tight;

laprint(1,'normconsist','options','factory');

