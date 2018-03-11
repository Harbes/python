%=========================================================================
%
%   Program to demonstrate two aspects of consistency
%
%=========================================================================

clear all
clc

%RandStream.setDefaultStream( RandStream('mt19937ar','seed',66) )

mu   = 0:0.1:20;
sig  = 4.0;

lnl1 = zeros( 0,length(mu));
lnl2 = zeros( 0,length(mu));
lnl3 = zeros( 0,length(mu));


% Population parameters 
mu_0  = 10.0;
sig_0 =  4.0;


% Sample T=5 
t = 5;
y = mu_0 + sig_0*randn(t,1);

disp( ['Sample mean (T=5)   =' num2str(mean(y)) ]);


for i = 1:length(mu)
    
    lnl1(i) = mean(-0.5*log(2*pi*sig^2) - 0.5*(y-mu(i)).^2./sig^2);

end

% Sample T=20
t = 20;
y = mu_0 + sig_0*randn(t,1);

disp( ['Sample mean (T=20)   =' num2str(mean(y)) ]);

for i = 1:length(mu)
    
    lnl2(i) = mean(-0.5*log(2*pi*sig^2) - 0.5*(y-mu(i)).^2./sig^2);

end




% Sample T=500
t = 500;
y = mu_0 + sig_0*randn(t,1);

disp( ['Sample mean (T=500)   =' num2str(mean(y)) ]);

for i = 1:length(mu)
    
    lnl3(i) = mean(-0.5*log(2*pi*sig^2) - 0.5*(y-mu(i)).^2./sig^2);

end



% Compute population log-likelihood  

e_lnl = -0.5*log(2*pi*sig_0^2) - 0.5 - 0.5*(mu - mu_0).^2/sig_0^2;

tmp = -0.5*(log(2*pi*sig_0^2) + 1);
disp( ['Maximum value at theta = theta_0 =' num2str(tmp) ] );


%********************************************************************
%***
%***     Generate graph
%***
%********************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

plot(mu,e_lnl,'-k',mu,lnl1,':k',mu,lnl2,'-.k',mu,lnl3,'--k','LineWidth',0.75);
axis([2 15 -3.5 -2.5])
ylabel('$A(\theta)$')
xlabel('$\mu$')

% Print TEX file
laprint(1,'twotypesconsistency','options','factory');

