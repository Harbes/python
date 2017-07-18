
%=========================================================================
%
%     Program to demonstrate the consistency property of MLE for the
%     location parameter (theta) of the Cauchy distribution.
%
%     For this example the median is the MLE and is thus a consistent estimator
%     but the sample mean is an inconsistent estimator.
%
%=========================================================================
clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234) )


theta  = 1;                             % Loaction parameter 
nu     = 1;                             % Cauchy = Student t nu =1 dof
t   = 500;

% Generate sample means and medians from sample of size t=1,2,3,...,t        

ybar = zeros( t,1 );
ymed = zeros( t,1 );

for i = 1:t

    y       = theta + trnd(nu,[i 1]);
    ybar(i) = mean(y);
    ymed(i) = median(y);
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

subplot(1,2,1)
plot(1:1:t,theta.*ones(t,1),'-k',1:1:t,ybar,'-k');
title('(a) Mean');
ylabel('$\widehat{\theta}$');
xlabel('Progressive Sample Size');
box off;
axis tight;

subplot(1,2,2)
plot(1:1:t,theta.*ones(t,1),'-k',1:1:t,ymed,'-k');
title('(b) Median');
ylabel('$\widehat{\theta}$');
xlabel('Progressive Sample Size');
box off;
axis tight;




%laprint(1,'cauchyconsist','options','factory');

