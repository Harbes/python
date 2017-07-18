%========================================================================
%
%   Program to compute the White and Newey-West standard errors 
%   for the OLS estimator
%
%========================================================================

clear all;
clc;

% Data to reproduce the numbers in the chapter 
y = [ 3,2,4,6,10 ]';
x = [ 1,2,3,4,5 ]';

% Data for part (d) of exercise
%y = [ 2,5,8,9,11 ];       
%x = [ 2,4,6,5,8 ];

t = length(y);

% Estimate by OLS
x  = [ones(t,1)  x];
b  = x\y;
u  = y - x*b;
s2 = mean( u.^2 );

% Compute gradients at each t
g = bsxfun(@times,u,x/s2);



% Compute covariance based on Hessian
h = -x'*x/s2;

disp('Covariance matrix (Hessian) ');
disp(-h); 
disp('Standard errors   (Hessian) ');
disp(sqrt(diag(-h)));

% Compute OPG
 j = g'*g;

disp('Covariance matrix (OPG)');
disp(j); 
disp('Standard errors   (OPG)');
disp(sqrt(diag(j)));

% White covariance matrix
covw = inv(h)*j*inv(h);
 
disp('Covariance matrix (White)');
disp(covw);
disp('Standard errors   (White)');
disp(sqrt(diag(covw)));
disp(' ');

% Newey-West covariance matrix
p = floor( 4*(t/100)^(2/9) );

for i = 1:p
 
    gmat = g((i+1):t,:)'*g(1:(t-i),:);
    gmat = trimr(g,i,0)'*trimr(g,0,i);
    j    = j + (1.0 - i/(p+1))*(gmat + gmat');

end
covnw = inv(h)*j*inv(h);

disp('Covariance matrix (Newey-West)');
disp(covnw);
disp('Standard errors   (Newey-West)');
disp(sqrt(diag(covnw)));
disp(' ');


