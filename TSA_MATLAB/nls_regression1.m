%=========================================================================
%
%   Program to estimate a nonlinear regression model
%
%=========================================================================
function nls_regression1( )

    clear all
    clc
    
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) );

    t = 100;
    
    % Generate the data                   
    theta = 2;
    x = 2 + rand(t,1);
    u = randn(t,1);
    y = 1./( x - theta) + u;

    % Method of Scoring   
    start = 1.5;              
    g = mean( (y - 1./(x - start)) .* (1./(x - start).^2) );
    i = mean( (1./(x - start).^4) );

    theta1 = start + inv(i)*g;

    disp('Method of Scoring')
    disp(['Starting value of theta = ',num2str(start) ]);
    disp(['Updated value of theta  = ',num2str(theta1) ]);

	% Estimate the model using BGS and compute Hessian se   
    ops = optimset('LargeScale', 'off', 'Display', 'iter');     
    
    [bhat,~,~,~,~,hess] = fminunc(@(b) neglog(b,y,x),start,ops);
    
    vc = (1/t)*inv(hess);
    disp(' ');
    disp( ['BFGS estimate of theta  = ',num2str(bhat) ])
    disp( ['Std. error  of theta    = ',num2str(sqrt(vc)) ])


end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
% Log-likelihood function
%-------------------------------------------------------------------------
function lf = neglog(b,y,x)

    m   = 1./(x - b);
    s2  = 1;
    lt = - 0.5*log(2*pi) - 0.5*log(s2) - 0.5*((y - m).^2)/s2;
    lf = -mean(lt);

end


