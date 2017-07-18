%=========================================================================
%
%     Program to demonstrate the Lindberg-Levy central limit theorem
%     using the uniform distribution.
%
%=========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',12) )

t      = 10;
Ndraws = 5000;                            % Number of draws                                     
mu     = 1/2;                             % Population mean of the uniform distribution         
sig2   = 1/12;                            % Population variance of the uniform distribution     
u      = rand(t,Ndraws);                  % Generate uniform random numbers                 

% For the sample size of 2        
y     = u(1:2,: );                          
ybar1 = mean(y);
z1    = sqrt(2)*(ybar1 - mu)./sqrt(sig2);    % Standardized variable 

% For the sample size of 10          
y     = u(1:10,: );                           
ybar2 = mean(y);
z2    = sqrt(10)*(ybar2 - mu)./sqrt(sig2);    % Standardized variable 


%**************************************************************************
%**     
%**     Generate graph       
%**                                         
%**************************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;


subplot(2,2,1)
hist(z1,21);
title('(a) Distribution of $z$ (T = 2)');
ylabel('$f(z)$');
%xlabel('$z$');
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
box off;


subplot(2,2,2)
hist(z2,21);
title('(b) Distribution of $z$ (T = 10)');
ylabel('$f(z)$');
%xlabel('$z$');
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
box off;


subplot(2,2,3)
hist(ybar1,21);
title('(c) Distribution of $\bar{y}$ (T = 2)');
ylabel('$f(\bar{y})$');
%xlabel('$\bar{y}$');
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
box off;


subplot(2,2,4)
hist(ybar2,21);
title('(d) Distribution of $\bar{y}$ (T = 10)');
ylabel('$f(\bar{y})$');
%xlabel('$\bar{y}$');
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
box off;


% Print TEX file
laprint(1,'clt','options','factory');




