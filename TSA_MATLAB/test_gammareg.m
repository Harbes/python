%   Exercise 4.8

function ch4_8

rand('state',0)

%   Part (a) - simulation

beta0 = 1;
beta1 = 2;
T     = 2000;

x_t  = rand(T,1);
mu_t = beta0 + beta1.*x_t;
u_t  = rand(T,1);
y_t  = -mu_t.*log(1 - u_t);

%   Part (b) - unconstrained estimation

th0 = [1;1;1];
options = optimset('GradObj','off','Display','off');
[theta_1,logL_1] = fminsearch(@logl_u,th0,options,y_t,x_t);

fprintf('\n')
fprintf('Unconstrained estimation results \n')
fprintf('beta0_hat_1  = %3.4f \n',theta_1(1))
fprintf('beta1_hat_1  = %3.4f \n',theta_1(2))
fprintf('rho_hat_1    = %3.4f \n',theta_1(3))
fprintf('logL_1       = %3.4f \n',-logL_1)
fprintf('\n')

%   Part (c) - testing

%   (i)     Likelihood ratio test

th0 = [1;1;1];
options = optimset('GradObj','off','Display','off');
[theta_0,logL_0] = fminsearch(@logl_c,th0,options,y_t,x_t);
theta_0(3) = 1;

fprintf('\n')
fprintf('Constrained estimation results \n')
fprintf('beta0_hat_0  = %3.4f \n',theta_0(1))
fprintf('beta1_hat_0  = %3.4f \n',theta_0(2))
fprintf('rho_hat_1    = %3.4f \n',theta_0(3))
fprintf('logL_0       = %3.4f \n',-logL_0)
fprintf('\n')

LR = 2*(logL_0 - logL_1);
p  = 1 - chi2cdf(LR,1);

fprintf('Likelihood ratio test \n')
fprintf('LR stat      = %3.4f \n',LR)
fprintf('p-value      = %3.4f \n',p)
fprintf('\n')

%   (ii)    Wald test with Hessian computed from numerical derivatives

H_1 = -numhess(@logl_u,theta_1,y_t,x_t);

R = [0 0 1];
Q = 1;

fprintf('Hessian evaluated at theta_hat_1, determined numerically \n')
fprintf('H(th_hat_1)  = \n')
fprintf('               %4.4f  %4.4f %4.4f \n',H_1)
fprintf('\n')

WH = (R*theta_1 - Q)'*inv( R*inv(-H_1)*R' )*(R*theta_1 - Q);
pH = 1 - chi2cdf(WH,1);

fprintf('Wald tests with numerical derivatives \n')
fprintf('Using Hessian \n')
fprintf('Wald stat    = %3.4f \n',WH)
fprintf('p value      = %3.4f \n',pH)
fprintf('\n')

%   (iii)   LM test with OPG computed from numerical derivatives

Gt  = numgrad(@loglt_u,theta_0,y_t,x_t);
J_0 = Gt'*Gt;

fprintf('OPG evaluated at theta_hat_0, determined numerically \n')
fprintf('J(th_hat_0)  = \n')
fprintf('               %4.4f  %4.4f %4.4f \n',J_0)
fprintf('\n')

G_0 = sum(Gt)';

fprintf('gradient evaluated at theta_hat_0 \n')
fprintf('G(th_hat_0)  = \n')
fprintf('               %4.4f \n',G_0)
fprintf('\n')

LMJ      = G_0'*inv(J_0)*G_0;
pJ       = 1 - chi2cdf(LMJ,1);

fprintf('LM test with OPG matrix computed from numerical derivatives \n')
fprintf('LM stat      = %3.4f \n',LMJ)
fprintf('p value      = %3.4f \n',pJ)
fprintf('\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Sub-routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Define the unconstrained log-likelihood at each observation

function LogLt = loglt_u(theta,y_t,x_t)

beta0 = theta(1);
beta1 = theta(2);
rho   = theta(3);
mu_t  = beta0 + beta1*x_t;

LogLt = -log(gamma(rho)) - rho.*log(mu_t) + (rho-1).*log(y_t) -...
    y_t./mu_t;

%   -log likelihood summed over all observations

function LogL = logl_u(theta,y_t,x_t)

LogLt = loglt_u(theta,y_t,x_t);
LogL  = -sum(LogLt);

%   constrained log-likelihood at each observation

function LogLt = loglt_c(theta,y_t,x_t)

beta0 = theta(1);
beta1 = theta(2);
rho   = 1;
mu_t  = beta0 + beta1*x_t;

LogLt = -log(gamma(rho)) - rho.*log(mu_t) + (rho-1).*log(y_t) -...
    y_t./mu_t;

%   -log likelihood summed over all observations

function LogL = logl_c(theta,y_t,x_t)

LogLt = loglt_c(theta,y_t,x_t);
LogL  = -sum(LogLt);
