%=========================================================================
%
%     Monte Carlo program to demonstrate the asymptotical normality
%     of MLE with R = 5000 replications based on an exponential distribution
%     with theta = 1.
%
%=========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234) )

R     = 5000;                           
theta = 1;


%     For the sample size of T = 5      
%

t = 5;                             %   Sample size 
u = rand(t,R);                     %   Generate uniform random numbers         
y = -theta.*log(1 - u);            %   Generate realizations of y              

ybar = mean(y);
z1   = sqrt(t).*(ybar - theta)./sqrt(theta);  %   Standardized random variable    


%     For the sample size of T = 100      
%

t = 100;                            %  Sample size                     
u = rand(t,R);                      %  Generate uniform random numbers                 
y = -theta.*log(1 - u);             %  Generate realizations of y                      

ybar = mean(y);
z2   = sqrt(t)*(ybar - theta)./sqrt(theta);  %   Standardized random variable     


%**************************************************************************
%**
%**     Generate graph
%**
%**************************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

subplot(2,2,1:2)
s  = 0:0.001:12;
fy = exp(-s./theta)./theta;
 
plot(s,fy,'-k','LineWidth',0.75);
title('(a) Exponential distribution');
ylabel('$f(y)$');
xlabel('$y$');
ylim( [0, 1.5] )
xlim([0,6]);
set(gca,'XTick',0:2:6);
set(gca,'YTick',0:0.5:1.5);
box off;
%axis tight;



subplot(2,2,3)
hist(z1,21)
title('(b) $T=5$');
ylabel('$f(z)$');
xlabel('$z$');
xlim([-5,5]);
set(gca,'XTick',-4:1:4);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
box off;
%axis tight;


subplot(2,2,4)
hist(z2,21)
title('(c) $T=100$');
ylabel('$f(z)$');
xlabel('$z$');
xlim([-5,5]);
set(gca,'XTick',-4:1:4);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
box off;
%axis tight;
 

%laprint(1,'aympnorm','options','factory');