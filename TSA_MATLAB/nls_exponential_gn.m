%==========================================================================
%
%    Program to estimate an exponential model
%    using the GAUSS-NEWTON algorithm 
%     
%==========================================================================

clear all;
clc;       

RandStream.setDefaultStream( RandStream('mt19937ar','seed',123457) );


%  Simulate data
t   = 50;        			    
b0  = 1.0;
b1  = 0.05;
sig = 0.5;
u   = sig*randn(t,1);
x   = (1:1:t)';
y   = b0*exp( b1*x ) + u;

% GAUSS-NEWTON algorithm  

b = [ 0.1; 0.1];	 % Choose staring values       
crit  = 0.00001;     % Convergence criterion         
maxit = 20;          % Maximum number of iterations  

for i=1:maxit

    e  = y - b(1)*exp(b(2)*x);
    z1 = exp(b(2)*x);
    z2 = b(1)*x.*exp(b(2)*x);
    z  = [z1 z2];

    disp( ['Iteration        = ', num2str(i) ] );

    bchange = z\e;

    if (bchange'*bchange) < crit; 
        break; 
    else
        b = b + bchange;			 % Update new parameters   

	if i==maxit;
	    fprintf(1,'Failed to converge after iteration number %d\n', maxit);
	end

    end

end
s2    = (e'*e)/t;                       % Residual variance       
omega = s2*inv(z'*z);

disp( 'Parameter estimates, std errors and t-stats' );
disp( [b sqrt(diag(omega)) b./sqrt(diag(omega)) ] );

disp(' ');
disp('Estimated asymptotic covariance matrix') 
disp( omega );

