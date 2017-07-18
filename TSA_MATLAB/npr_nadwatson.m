%=========================================================================
%
%    Nonparametric kernel regression estimates of the numerical example 
%    
%=========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234) );

% Simulate the model
t  = 20;
ut = 0.1*randn(t,1);                % N(0,0.1^2) 
xt = -2 + 4*rand(t,1);              % U(-2,2) 
mx = 0.3*exp( -4*(xt + 1).^2 ) + 0.7*exp( -16*(xt - 1).^2 );  
yt  = mx + ut;

% Construct true conditional mean
x   = -2:0.01:2;
mxt = 0.3*exp( -4*(x + 1).^2 ) + 0.7*exp( -16*(x - 1).^2 );  

%  Compute kernel regression for h = 0.5  
h  = 0.2;
fx = zeros(length(x),1);
fxy = zeros(length(x),1);
for i = 1:length(x);
    z      = ((x(i) - xt)./h);    
    fx(i)  = mean( normpdf(z)./h );
    fxy(i) = mean( normpdf(z).*yt./h );
end
mx1 = fxy ./ fx;

% Compute kernel regression for h = 0.1     **/
h  = 0.02;
for i = 1:length(x);
    z      = ((x(i) - xt)./h);    
    fx(i)  = mean( normpdf(z)./h );
    fxy(i) = mean( normpdf(z).*yt./h );
end
mx2 = fxy ./ fx;


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
subplot(1,2,1);
plot(x,mx1,'-.k','LineWidth',1.0);
hold on;
plot(x,mxt,'-k','LineWidth',0.5);
plot(xt,yt,'.k');
hold off;
title('(a) Larger Bandwidth');
ylabel('$y_t, \, m(x_t)$');
xlabel('$x_t$');
box off;
hold off;

%--------------------------------------------------------%
% Panel (b)
subplot(1,2,2);
plot(x,mx2,'-.k','LineWidth',1.0);
hold on;
plot(x,mxt,'-k','LineWidth',0.5);
plot(xt,yt,'.k');
title('(b) Smaller Bandwidth');
ylabel('$y_t, \, m(x_t)$');
xlabel('$x_t$');
box off;

%laprint(1,'nadwatson','options','factory');

