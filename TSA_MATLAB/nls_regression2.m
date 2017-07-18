%=========================================================================
%
%   Program to estimate a nonlinear regression model
%
%=========================================================================
function nls_regression2( )

    clear all
    clc
    
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',12) )
    
    t = 100;
    
    % Parameters
    beta0 = 10;
    beta1 = 2;
    beta2 = 0.5;
    sig2  = 0.1;
    
    % Generate the data                   
    x = rand(t,1).^2;
    u = randn(t,1);
    y = (beta0 + beta1*x + sqrt(sig2)*u).^(1/beta2);

	% Estimate the model using BGS and compute Hessian se   
    start = [beta0 ; beta1; beta2 ];
    ops   = optimset('LargeScale', 'off', 'Display', 'iter');     
    
    [bhat,~,~,~,~,hess] = fminunc(@(b) neglog(b,y,x),start,ops);
    
    vc = (1/t)*inv(hess);
    disp(' ');
    disp( '   Estimates  ' )
    disp( bhat )

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
% Log-likelihood function
%-------------------------------------------------------------------------
function lf = neglog(b,y,x)

     t = length(y);
     m = b(1) + b(2)*x;
     u = y.^b(3) - m;

     % Concentrate the likelihood
     s2 = u'*u/t;             
     lt = - 0.5*log(2*pi) - 0.5*log(s2) + log(b(3)) + (b(3)-1)*log(y) - 0.5*(u.^2)/s2;
     lf = -mean(real(lt));

end

