%=========================================================================
%
%     Program to find the MLEs using grid search methods
%
%=========================================================================

clear all;
clc;

%RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) )


% Simulate the model  

x = [ 1; 2; 4; 5; 8; ];

beta = 1.0;
sig2 = 4.0;

t = size(x,1);
y = beta*x + sqrt(sig2)*randn(t,1);


% Grid search on gradient sig2 = 4      

sig2 = 4.0;
beta = (0.5:0.1:1.5)';
g1   = zeros(11,1);

for i=1:11
    g1(i) = mean( (y - beta(i)*x).*x )' / sig2;
end

figure(1)

subplot(1,2,1);
plot(beta,g1,beta,zeros(length(beta)),'-k');
title('Gradient: sig2 = 4.0');


% Grid search on gradient sig2 = 3.5   

sig2 = 3.5;
beta = (0.5:0.1:1.5)';
g2   = zeros(11,1);

for i=1:11
    g2(i) = mean( (y - beta(i)*x).*x )' / sig2;
end

subplot(1,2,2)
plot(beta,g2,beta,zeros(length(beta)),'-k');
title('Gradient: sig2 = 3.5');



