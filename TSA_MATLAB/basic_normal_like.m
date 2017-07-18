%=========================================================================
%
%   Program to estimate the parameters of a normal distribution by
%   maximum likelihood and plot the log of the likelihood function.
%
%=========================================================================

clear all;
clc;

% Read in data
y = [ 5, -1, 3, 0, 2, 3]';
t = length(y);       	% Define the sample size  

% Compute the MLEs

mu_mle   = mean(y);
sig2_mle = mean( (y - mu_mle).^2 );

disp('MLE(mu)');
disp('-------');
disp(mu_mle);

disp('MLE(sig2)');
disp('-------');
disp(sig2_mle);

%**************************************************************************
%**
%**     Generate graph
%**
%**************************************************************************

mu   =  1.0:0.05:3.00;
sig2 =  3.0:0.05:5.0;

lnl = zeros( length(mu),length(sig2) );
for i = 1:length(mu);
    for j = 1:length(sig2)
        
        lnl(i,j) = -0.5*t*log(2*pi)-0.5*t*log(sig2(j))-sum( (y - mu(i)).^2/sig2(j) );
    end
end
        
% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

mesh(mu,sig2,lnl,'EdgeColor','black')

xlabel('$\mu$');
ylabel('$\sigma^2$');
zlabel('$\mathrm{ln} L( \mu,\, \sigma^2)$');
set(gca,'ztick',[])
axis tight
grid 'off'
box 'off'

%saveas(gcf,'D:\Stan\Current\Vancebook\mle-basic\Graphics\m_c1fig5.eps','eps')
laprint(1,'normlike','options','factory');

