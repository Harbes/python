%=========================================================================
%
%     Program to estimate a nonlinear consumption function
%     using the GAUSS-NEWTON algorithm
%
%     U.S. data on real consumption and real disposable income(2005 $)
%     1960:Q1 to 2009:Q4 200 observations
%=========================================================================

clear all;
clc;

% Load data
load USdata.mat

yt = inc;
ct = cons;
t = length( yt );

% Estimate the linear model to use as initial starting values  
b = [ones(t,1) yt]\ct;

% GAUSS-NEWTON algorithm   
alpha = b(1);
beta  = b(2);
gam   = 1.00;

crit  = 0.00001;                        % Convergence criterion         
maxit = 20;                             % Maximum number of iterations  

for i=1:maxit

    et  = ct - alpha - beta*yt.^gam;           
    z1t = ones(t,1);                         
    z2t = yt.^gam;
    z3t = beta*log(yt).*yt.^gam;

    disp( ' ' )
    disp(['Iteration        = ' num2str(i)]) ;
    disp(['Parameters       = ' num2str([alpha beta gam]) ] );
    disp('---------------------------------------------------------------');
    y       = et;
    x       = [z1t z2t z3t];
    bchange = x\y;

    if (bchange'*bchange) < crit 
    	break; 
    else
        alpha = alpha + bchange(1);     % Update new parameters   
        beta  = beta  + bchange(2);
        gam   = gam   + bchange(3);

    	if i == maxit
	    disp(['Failed to converge after iteration: ', num2str(maxit)] );
    	end
    end
end


ssr  = y'*y;
sig2 = ssr/t;
disp( ['Sum squared residuals = ' num2str(ssr)] ); 
disp( ['Residual variance     = ' num2str(sig2) ] );


disp( 'Information matrix' );
disp('---------------------');
z = [z1t z2t z3t];
im = z'*z/sig2;
disp( im );

disp('Estimated asymptotic covariance matrix:\n');
disp( inv(im) );

bhat = [alpha beta gam]';
se   = sqrt(diag(inv(im)));
disp('Parameter estimates, standard errors and t-stats');
disp( [bhat se bhat./se]);


