%=========================================================================
%
%    Program to demonstrate the convergence properties of the average
%    log likelihood to its expectation assuming a normal distribution.
%
%=========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',1) )

t    = 1000;    
mu   = 1;               % Population mean
sig2 = 2;               % Population variance  

% Generate data from a normal distribution        
y = mu + sqrt(sig2)*randn(t,1);

% Compute progressive sample means of the log-likelihood         

lnf = -0.5*log(2*pi) - 0.5*log(sig2) - 0.5*(y - mu).^2/sig2;
tt  = (1:1:t)';

s_lnf = cumsum(lnf)./tt;
e_lnf = -0.5*log(2*pi) - 0.5*log(sig2) - 0.5;

disp('Population mean');
disp('----------------');
disp(e_lnf);

disp('Sample mean');
disp('-----------');
disp(mean(s_lnf));


%**************************************************************************
%**     
%**     Generate graph       
%**                                         
%**************************************************************************

%Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

plot(tt,e_lnf.*ones(t,1),':k','LineWidth',1.5);
hold on
plot(tt,s_lnf,'-k','LineWidth',0.75);

ylabel('$A(\theta)$');
xlabel('T');

ylim([-1.9 -1.6])
set(gca,'YTick',[-1.8 -1.7 -1.6])
box off;

% Print TEX file
laprint(1,'avelogL','options','factory');

