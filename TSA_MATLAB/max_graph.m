
%=========================================================================
%
%     Program to find the MLEs using graphical methods
%
%=========================================================================

clear all;
clc;

%RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) )


% Simulate the model   

x    = [ 1; 2; 4; 5; 8 ];
beta = 1.0;
sig2 = 4.0;
t    = length( x );
y    = beta*x + sqrt(sig2)*randn(t,1);


% Plot log-likelihood sig2 = 4

sig2 = 4.0;
beta = (0.0:0.1:2.0)';
a1   = zeros(21,1);

for i=1:length(beta)
    z = ( y - beta(i)*x )/sqrt(sig2);
    a1(i) = -0.5*log(2*pi) - 0.5*log(sig2) - 0.5*mean( z.^2 );
end
figure;
subplot(2,2,1);
plot(beta,a1);
title('Log-likelihood: sig2 = 4.0');


% Plot log-likelihood sig2 = 3.5

sig2 = 3.5;
beta = (0.0:0.1:2.0)';
a2   = zeros(21,1);

for i=1:21
    z = ( y - beta(i)*x)/sqrt(sig2);
    a2(i) = -0.5*log(2*pi) - 0.5*log(sig2) - 0.5*mean( z.^2 );
end

subplot(2,2,2);
plot(beta,a2);
title('Log-likelihood: sig2 = 3.5');


% Plot log-likelihood beta = 0.9

beta = 1.0;
sig2 = (1.0:0.5:11)';
a3   = zeros(21,1);

for i=1:21
    z = ( y - beta*x)/sqrt(sig2(i));
    a3(i) = -0.5*log(2*pi) - 0.5*log(sig2(i)) - 0.5*mean( z.^2 );
end

subplot(2,2,3);
plot(sig2,a3);
title('Log-likelihood: beta = 1.0');
axis tight

% Plot log-likelihood beta = 0.9

beta = 0.9;
sig2 = (1.0:0.5:11)';
a4   = zeros(21,1);

for i=1:21
    z = ( y - beta*x)/sqrt(sig2(i));
    a4(i) = -0.5*log(2*pi) - 0.5*log(sig2(i)) - 0.5*mean( z.^2 );
end

subplot(2,2,4);
plot(sig2,a4);
title('Log-likelihood: beta = 0.9');
axis tight

