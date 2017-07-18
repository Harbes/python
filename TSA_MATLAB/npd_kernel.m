% ========================================================================
%
%      Gaussian kernel  
%
% ========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',123457) );

xi = [ 2, 3, 5, 7, 8, 8, 9, 9, 10, 10, 11, 11, 11, 12, 12, 14, 15, 17, 18, 20]';
n  = size(xi,1);

% Compute kernel with bandwith h = 1       
h  = 1.0;
x  = (5:5:15)';
w1 = (normpdf(repmat(x,1,n)-repmat(xi',3,1))/h)';
f1 = sum(w1)/(n*h); 

% Compute kernel with bandwith h = 2       
h  = 2.0;
w2 = (normpdf(repmat(x,1,n)-repmat(xi',3,1))/h)';
f2 = sum(w2)/(n*h); 

disp( [w1 w2]);

% Plot the results
plot(x,f1,'-r',x,f2,'-g');

