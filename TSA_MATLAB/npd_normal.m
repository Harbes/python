% ========================================================================
%
%     Normal distribution example using different bandwidths.
%
% ========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',1) );

% Simulate the data from a normal distribution      
t   = 201;       
mu  = 0;
sig = 3;
yt  = mu + sig*randn(t,1);

% Generate the population (normal) distribution      
y = -10:0.1:10;
f_norm = normpdf( ( y - mu )/sig )/sig;

% Estimate the nonparametric density with h = 1      
f1 = zeros(length(y),1);
h  = 1;
for i = 1:length(y);
    z     = (y(i) - yt)./h;    
    f1(i) = mean( normpdf(z)./h );
end

% Estimate the nonparametric density with h = 2       
f2 = zeros(length(y),1);
h  = 2;
for i = 1:length(y);
    z     = (y(i) - yt)./h;    
    f2(i) = mean( normpdf(z)./h );
end

% Estimate the nonparametric density with h based on rot  
f3 = zeros(length(y),1);
h  = 1.06*std(yt)*t^(-1/5);
for i = 1:length(y);
    z     = (y(i) - yt)./h;    
    f3(i) = mean( normpdf(z)./h );
end

% Estimate the parametric density assuming normality     
m = mean(yt);
s = std(yt);
f_para = normpdf( ( y - m )/s )/s;

%**************************************************************************
%**
%**     Generate graph
%**
%**************************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

%--------------------------------------------------------%
% Panel (a)
subplot(1,3,1);
plot(y,f1,'-.r','LineWidth',1.0);
hold on;
plot(y,f_norm','-k','LineWidth',0.5);
plot(y,f_para','-k','LineWidth',0.2);
hold off;
title('(a)');
ylabel('$f(y)$');
xlabel('$y$');
box off;

%--------------------------------------------------------%
% Panel (b)
subplot(1,3,2);
plot(y,f2,'-.r','LineWidth',1.0);
hold on;
plot(y,f_norm','-k','LineWidth',0.5);
plot(y,f_para','-k','LineWidth',0.2);
hold off;
title('(b) h=2 ');
ylabel('$f(y_t)$');
xlabel('$y$');

box off;


%--------------------------------------------------------%
% Panel (c)
subplot(1,3,3);
plot(y,f3,'-.r','LineWidth',1.0);
hold on;
plot(y,f_norm','-k','LineWidth',0.5);
plot(y,f_para','-k','LineWidth',0.2);
hold off;
title('(c) h based on rot ');
ylabel('$f(y_t)$');
xlabel('$y$');

box off;
