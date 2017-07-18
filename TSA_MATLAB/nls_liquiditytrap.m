%=========================================================================
%
%   Program to test a liquidity trap for the United States
%
%=========================================================================
function nls_liquiditytrap( )

    clear all
    clc

    % Load data for the United States: January 1959 to December 2011   
    % Variables are
    %       m2
    %       gdp (real)
    %       cpi
    %       interest
    load us_liquiditytrap

    % Construct variables 
    y  = log(m2./cpi);        
    x1 = log(gdp);           
    x2 = interest/100;       

    t = length(y);

    % Estimate linear model by OLS         
    x    = [ones(t,1)    x1   x2 ];
    bols = x\y;
    u    = y - x*bols;
    s2 = u'*u/t;

    % LM test (2-step) 
    % Intercept is used in the construction of the test 
    % so the restricted model has an intercept. **/

    % Stage 1 regression
    x = [ones(t,1)   x1   1./x2 ];         
    b1 = x\y;                
    u = y - x*b1;      

    % Stage 2 regression
    z  = [x   1./(x2.^2) ];  
    b2 = z\u;
    v  = u - z*b2;                      
    r2 = 1 - (v'*v)/(u'*u);                   
    lm = t*r2;              

    disp(' ')
    disp(['LM test (2-step) with intercept  = ',num2str(lm) ]);
    disp(['p-value                          = ',num2str(1-chi2cdf(lm,1)) ]); 
    
    % Estimate model by MLE
    start = bols([2 1 3]);    % Note change of order of OLS estimates
    ops   = optimset('LargeScale','off','Display','off');
    
    [ bhat,~,~,~,~,hess ] = fminunc(@(b) neglog(b,y,x1,x2),start);

    % fminunc hess doesn't seem correct
    hess = numhess(@neglog,bhat,y,x1,x2 );

    vc = (1/t)*inv(hess);
    
    % Wald test
    wd = (bhat(3) - 0)^2/vc(3,3);

    disp(' ')
    disp(['Wald test        = ',num2str(wd) ]);
    disp(['p-value          = ',num2str(1-chi2cdf(wd,1)) ]);
 
    
end
%
%--------------------------- Functions -----------------------------------
% 
%------------------------------------------------------------------------- 
%   Negative log-likelihood function
%------------------------------------------------------------------------- 
function lf = neglog(b,y,x1,x2)

        t   = length(y);
        m   = b(1)*x1 +  b(2)./(x2 -  b(3));
        u   = y - m;
        s2  = u'*u/t;
        lnl = - 0.5*log(2*pi) - 0.5*log(s2) - 0.5*((y - m).^2)/s2;

        lf = -mean( lnl );
end
