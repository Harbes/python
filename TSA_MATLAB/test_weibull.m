%=========================================================================
%
%     Program to do LR, Wald and LM test of Weibull Distribution
%
%=========================================================================
function test_weibull( )

    clear all
    clc

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',12) )
    
    % Simulate data
    t     = 20;
    alpha = 1;
    beta  = 2;

    %y = wblrnd(alpha,beta,t,1);
    
    % or load data
    y = [0.293, 0.589, 1.374, 0.954, 0.608, 1.199, 1.464, ...
        0.383, 1.743, 0.022, 0.719, 0.949, 1.888, 0.754, ...
        0.873, 0.515, 1.049, 1.506, 1.090, 1.644]';

    
    bstart = [alpha; beta];
    options = optimset('Display','off',...
                       'LargeScale','off');

    % Unconstrained estimation
    [bu,logL1] = fminunc(@(b) loglu(b,y),bstart,options);
    logL1 = -logL1;
    
    disp(' ');
    disp('Unconstrained estimation results');
    disp( ['alpha  = ' num2str(bu(1)) ] );
    disp( ['beta   = ' num2str(bu(2)) ] );
    disp( ['log L  = ' num2str(t*logL1) ] );
    fprintf('\n')

    % Constrained estimation
    [b0,logL0] = fminunc(@(b) loglc(b,y),1,options);
    logL0 = -logL0;
    
    disp(' ')
    disp('Constrained estimation results');
    disp( ['alpha  = ' num2str(b0) ] );
    disp( ['beta   = ' num2str(1) ] );
    disp( ['log L  = ' num2str(t*logL0) ] );

    % Likelihood ratio test
    lr = -2*t*(logL0 - logL1);
    p  = 1 - chi2cdf(lr,1);

    disp('Likelihood ratio test')
    disp( ['LR stat      = ' num2str(lr) ] );
    disp( ['p-value      = ' num2str(p) ] );
    disp( ' ' );
 
    % Wald test
    h = numhess(@loglu,bu,y);
    
    disp('Hessian evaluated at unconstrained estimates');
    disp( h );
    disp( ' ' );

    r = [ 0 1 ];
    q = 1;
    w = t*(r*bu - q)'*inv(r*h*r')*(r*bu - q);
    p = 1 - chi2cdf(w,1);
    
    disp('Wald test')
    disp( ['Wald stat    = ' num2str(w) ] );
    disp( ['p-value      = ' num2str(p) ] );
    disp( ' ' );

    % LM test
    th0  = [b0; 1];
    gmat = numgrad(@logltu,th0,y);
	g    = mean(gmat)';
    j    = gmat'*gmat/t;
    lm   = t*g'*inv(j)*g;
    p    = 1 - chi2cdf(lm,1);

    
    disp('LM test')
    disp( ['LM stat    = ' num2str(lm) ] );
    disp( ['p-value    = ' num2str(p) ] );

end

%
%--------------------------- Functions -----------------------------------
% 

%-------------------------------------------------------------------------
%  Unconstrained likelihood function 
%-------------------------------------------------------------------------
function logl = loglu(b,y)

    logl = -mean( logltu(b,y) );

end

%-------------------------------------------------------------------------
%  Unconstrained likelihood function at each observation
%-------------------------------------------------------------------------
function lt = logltu(b,y)

    alpha = b(1);
    beta  = b(2);
    lt = log(alpha)+log(beta)+(beta-1)*log(y)- alpha*y.^(beta);
    
end

%-------------------------------------------------------------------------
%  Constrained likelihood function
%-------------------------------------------------------------------------
function logl = loglc(b,y)

    logl = -mean(logltc(b,y));

end

%-------------------------------------------------------------------------
%  Constrained likelihood function at each observation
%-------------------------------------------------------------------------
function lt = logltc(alpha,y)

    beta = 1;
    lt = log(alpha)+log(beta)+(beta-1)*log(y)- alpha*y.^(beta);


end
