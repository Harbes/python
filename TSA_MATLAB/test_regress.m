%=========================================================================
%
%     Program to do LR, Wald and LM tests of a simple regression model
%
%=========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',12) )

t    = 5;
x    = [1,2,4,5,8]';
beta = 1;
sig  = 2;

% Generate data
u = sig*randn(t,1);
y = x*beta + u;

% Constrained estimation
bhat0 = 0.0;
sig20 = mean((y - x*bhat0).^2);

logL0  = - 1/2*log(2*pi) - 1/2*log(sig20)...
               - 1/(2*sig20*t)*sum((y - x*bhat0).^2);

fprintf('Constrained estimation results \n')
fprintf('beta_hat_0   = %3.4f \n',bhat0)
fprintf('sigma2_hat_0 = %3.4f \n',sig20)
fprintf('logL_0       = %3.4f \n',logL0)
fprintf('\n')
           
% Unconstrained estimation

bhat1 = x\y;
err2  = (y - x*bhat1).^2;
sig21 = mean(err2);
 

logL1 = - 1/2*log(2*pi) - 1/2*log(sig21)...
               - 1/(2*sig21*t)*sum((y - x*bhat1).^2);

fprintf('Unconstrained estimation results \n')
fprintf('beta_hat_1   = %3.4f \n',bhat1)
fprintf('sigma2_hat_1 = %3.4f \n',sig21)
fprintf('logL_1       = %3.4f \n',logL1)
fprintf('\n')

% Likelihood ratio test
lr  = -2*t*(logL0 - logL1);
pv  = 1-chi2cdf(lr,1);

disp('Likelihood ratio test')
disp( ['LR stat      =  ' num2str(lr) ] );
disp( ['p-value      =  ' num2str(pv) ] );
disp(' ')

% Wald test
r      = [ 1  0 ];
q      =  0 ;
theta1 = [bhat1; sig21];
%omega1 = sig21*inv(x'*x);
i = [ (x'*x)/(t*sig21)  0; 0 1/(2*(sig21)^2) ];
w      = t*(r*theta1 - q)'*inv(r*i*r')*(r*theta1 - q);
pv     = chi2cdf(w,1);

disp('Wald test')
disp( ['Unconstrained estimates = ' num2str( [bhat1 sig21 ] ) ] );
disp( ['Wald stat               =  ' num2str(w) ] );
disp( ['p-value                 =  ' num2str(pv) ] );
disp(' ')

% LM test analytic derivatives
g = [ x'*(y-x*bhat0)/(t*sig20); ...
      -(1/(2*sig20))+(1/(2*(sig20)^2*t))*sum( (y - x*bhat0).^2 ) ];
i = [ (x'*x)/(t*sig20)  0; 0 1/(2*(sig20)^2) ];
    
lm = t*g'*inv(i)*g;  
pv = chi2cdf(lm,1);


disp('LM test')
disp( ['Constrained estimates   = ' num2str( [0 sig20 ] ) ] );
disp( ['LM stat                 =  ' num2str(lm) ] );
disp( ['p-value                 =  ' num2str(pv) ] );
disp(' ')


      