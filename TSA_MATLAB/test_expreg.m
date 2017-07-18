%   Exercise 4.7

function exp_reg()

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

th0 = [1; 1];
options = optimset('GradObj','on','Display','off');
[theta_1,logL_1] = fminsearch(@logl_u,th0,options,y_t,x_t);

fprintf('\n')
fprintf('Unconstrained estimation results \n')
fprintf('beta0_hat_1  = %3.4f \n',theta_1(1))
fprintf('beta1_hat_1  = %3.4f \n',theta_1(2))
fprintf('logL_1       = %3.4f \n',-logL_1)
fprintf('\n')

%   Part (c) - constrained estimation

th0 = [1; 1];
options = optimset('GradObj','off','Display','off');
[theta_0,logL_0] = fminsearch(@logl_c,th0,options,y_t,x_t);
theta_0(2) = 0;

fprintf('\n')
fprintf('Constrained estimation results \n')
fprintf('beta0_hat_0  = %3.4f \n',theta_0(1))
fprintf('beta1_hat_0  = %3.4f \n',theta_0(2))
fprintf('logL_0       = %3.4f \n',-logL_0)
fprintf('\n')

%   Part (d) - tests

%   (i)     Likelihood ratio test

LR = 2*(logL_0 - logL_1);
p  = 1 - chi2cdf(LR,1);

fprintf('Likelihood ratio test \n')
fprintf('LR stat      = %3.4f \n',LR)
fprintf('p-value      = %3.4f \n',p)
fprintf('\n')

%   (ii)	Wald test with analytic derivatives

H_1 = HExp(theta_1,y_t,x_t);
I_1 = IExp(theta_1,y_t,x_t);
Gt  = GtExp(theta_1,y_t,x_t);
J_1 = Gt'*Gt;

fprintf('Analytically determined matrices for covariance \n')
fprintf('Hessian evaluated at theta_hat_1 \n')
fprintf('H(th_hat_1)  = \n')
fprintf('               %3.4f  %3.4f \n',H_1)
fprintf('Information matrix evaluated at theta_hat_1 \n')
fprintf('I(th_hat_1)  = \n')
fprintf('               %3.4f  %3.4f \n',I_1)
fprintf('OPG evaluated at theta_hat_1 \n')
fprintf('J(th_hat_1)  = \n')
fprintf('               %3.4f  %3.4f \n',J_1)
fprintf('\n')

R = [0 1];
Q = 0;

WH = (R*theta_1 - Q)'*inv( R*inv(-H_1)*R' )*(R*theta_1 - Q);
WI = (R*theta_1 - Q)'*inv( R*inv( I_1)*R' )*(R*theta_1 - Q);
WJ = (R*theta_1 - Q)'*inv( R*inv( J_1)*R' )*(R*theta_1 - Q);

pH = 1 - chi2cdf(WH,1);
pI = 1 - chi2cdf(WI,1);
pJ = 1 - chi2cdf(WJ,1);

fprintf('Wald tests with analytical derivatives \n')
fprintf('Using Hessian \n')
fprintf('Wald stat    = %3.4f \n',WH)
fprintf('p value      = %3.4f \n',pH)
fprintf('Using information matrix \n')
fprintf('Wald stat    = %3.4f \n',WI)
fprintf('p value      = %3.4f \n',pI)
fprintf('Using OPG \n')
fprintf('Wald stat    = %3.4f \n',WJ)
fprintf('p value      = %3.4f \n',pJ)
fprintf('\n')

%   (ii)	Wald test with numerical derivatives

Gt  = numgrad(@loglt_u,theta_1,y_t,x_t);
J_1 = Gt'*Gt;
H_1 = -numhess(@logl_u,theta_1,y_t,x_t);


fprintf('Numerically determined matrices for covariance \n')
fprintf('Hessian evaluated at theta_hat_1 \n')
fprintf('H(th_hat_1)  = \n')
fprintf('               %3.4f  %3.4f \n',H_1)
fprintf('OPG evaluated at theta_hat_1 \n')
fprintf('J(th_hat_1)  = \n')
fprintf('               %3.4f  %3.4f \n',J_1)
fprintf('\n')

WH = (R*theta_1 - Q)'*inv( R*inv(-H_1)*R' )*(R*theta_1 - Q);
WJ = (R*theta_1 - Q)'*inv( R*inv( J_1)*R' )*(R*theta_1 - Q);

pH = 1 - chi2cdf(WH,1);
pJ = 1 - chi2cdf(WJ,1);

fprintf('Wald tests with numerical derivatives \n')
fprintf('Using Hessian \n')
fprintf('Wald stat    = %3.4f \n',WH)
fprintf('p value      = %3.4f \n',pH)
fprintf('Using OPG \n')
fprintf('Wald stat    = %3.4f \n',WJ)
fprintf('p value      = %3.4f \n',pJ)
fprintf('\n')

%   (iii)	LR tests with analytic derivatives

H_0 = HExp(theta_0,y_t,x_t);
I_0 = IExp(theta_0,y_t,x_t);
Gt  = GtExp(theta_0,y_t,x_t);
J_0 = Gt'*Gt;
G_0 = GExp(theta_0,y_t,x_t)';

fprintf('Analytically determined results \n')
fprintf('gradient evaluated at theta_hat_0 \n')
fprintf('G(th_hat_0)  = \n')
fprintf('               %3.4f \n',G_0)
fprintf('\n')

fprintf('Hessian evaluated at theta_hat_0 \n')
fprintf('H(th_hat_0)  = \n')
fprintf('               %3.4f  %3.4f \n',H_0)
fprintf('Information matrix evaluated at theta_hat_0 \n')
fprintf('I(th_hat_0)  = \n')
fprintf('               %3.4f  %3.4f \n',I_0)
fprintf('OPG evaluated at theta_hat_0 \n')
fprintf('J(th_hat_0)  = \n')
fprintf('               %3.4f  %3.4f \n',J_0)
fprintf('\n')

LMH      = G_0'*inv(-H_0)*G_0;
LMI      = G_0'*inv(I_0)*G_0;
LMJ      = G_0'*inv(J_0)*G_0;

pH       = 1 - chi2cdf(LMH,1);
pI       = 1 - chi2cdf(LMI,1);
pJ       = 1 - chi2cdf(LMJ,1);

fprintf('LM tests with analytical derivatives \n')
fprintf('Using Hessian \n')
fprintf('LM stat      = %3.4f \n',LMH)
fprintf('p value      = %3.4f \n',pH)
fprintf('Using information matrix \n')
fprintf('LM stat      = %3.4f \n',LMI)
fprintf('p value      = %3.4f \n',pI)
fprintf('Using OPG \n')
fprintf('LM stat      = %3.4f \n',LMJ)
fprintf('p value      = %3.4f \n',pJ)
fprintf('\n')

%   (iii)   LR tests with numerical derivatives

Gt  = numgrad(@loglt_u,theta_0,y_t,x_t);
J_0 = Gt'*Gt;
H_0 = -numhess(@logl_u,theta_0,y_t,x_t);
G_0 = sum(Gt)';

fprintf('Numerically determined results \n')
fprintf('gradient evaluated at theta_hat_0 \n')
fprintf('G(th_hat_0)  = \n')
fprintf('               %3.4f \n',G_0)
fprintf('\n')
fprintf('Hessian evaluated at theta_hat_0 \n')
fprintf('H(th_hat_0)  = \n')
fprintf('               %3.4f  %3.4f \n',H_0)
fprintf('OPG evaluated at theta_hat_0 \n')
fprintf('J(th_hat_0)  = \n')
fprintf('               %3.4f  %3.4f \n',J_0)
fprintf('\n')

LMH      = G_0'*inv(-H_0)*G_0;
LMJ      = G_0'*inv(J_0)*G_0;

pH       = 1 - chi2cdf(LMH,1);
pJ       = 1 - chi2cdf(LMJ,1);

fprintf('LM tests with numerical derivatives \n')
fprintf('Using Hessian \n')
fprintf('LM stat      = %3.4f \n',LMH)
fprintf('p value      = %3.4f \n',pH)
fprintf('Using OPG \n')
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
mu_t  = beta0 + beta1*x_t;

LogLt = -log(mu_t) - y_t./mu_t;

%   -log likelihood summed over all observations

function LogL = logl_u(theta,y_t,x_t)

LogLt = loglt_u(theta,y_t,x_t);
LogL  = -sum(LogLt);

%   analytic gradient at each observation

function Gt = GtExp(theta,y_t,x_t)

T = size(y_t,1);

beta0 = theta(1);
beta1 = theta(2);
mu_t  = beta0 + beta1*x_t;

Gt      = zeros(T,2);
Gt(:,1) = -mu_t.^(-1) + y_t.*mu_t.^(-2);
Gt(:,2) = x_t.*Gt(:,1);

%   analytic gradient summed over all observations

function G = GExp(theta,y_t,x_t)

Gt = GtExp(theta,y_t,x_t);
G  = sum(Gt);

%   analytic hessian

function H = HExp(theta,y_t,x_t)

beta0 = theta(1);
beta1 = theta(2);
mu_t  = beta0 + beta1*x_t;

H      = zeros(2,2);
H(1,1) = sum(mu_t.^(-2) - 2.*y_t.*mu_t.^(-3));
H(1,2) = sum( (mu_t.^(-2) - 2.*y_t.*mu_t.^(-3)).*x_t );
H(2,1) = H(1,2);
H(2,2) = sum( (mu_t.^(-2) - 2.*y_t.*mu_t.^(-3)).*x_t.^2 );

%   analytic information matrix

function I = IExp(theta,y_t,x_t)

beta0 = theta(1);
beta1 = theta(2);
mu_t  = beta0 + beta1*x_t;

I      = zeros(2,2);
I(1,1) = sum(mu_t.^(-2));
I(1,2) = sum( (mu_t.^(-2)).*x_t );
I(2,1) = I(1,2);
I(2,2) = sum( (mu_t.^(-2)).*x_t.^2 );

%   constrained log-likelihood at each observation

function LogLt = loglt_c(theta,y_t,x_t)

beta0 = theta(1);
beta1 = 0;
mu_t  = beta0 + beta1*x_t;

LogLt = -log(mu_t) - y_t./mu_t;

%   -log likelihood summed over all observations

function LogL = logl_c(theta,y_t,x_t)

LogLt = loglt_c(theta,y_t,x_t);
LogL  = -sum(LogLt);

