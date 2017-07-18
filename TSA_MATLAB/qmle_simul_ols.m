%========================================================================
%
%      Simulation of data from an autoregressive model with heteroskedasticity
%      Estimate by OLS and generate various standard errors
%
%========================================================================

function qmle_simul_ols(  )

    clear all
    clc
   
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',12) );

    t     = 500;                
    beta1 = 1.0;   
    beta2 = 0.5; 
    rho   = 0.5;
    gam1  = 1.0;  
    gam2  = 0.5;

    % Construct exogenous variables 
    tmp = (1:1:t)';
    x = 0.01*tmp + 10*rand(t,1)*cos( 2*pi*t/40 );        

    % Simulate the AR(1) regression model with heteroskedasticity    

    y = randn(t,1);
    u = randn(t,1);
    v = randn(t,1);

    for i = 2:t

        sig2 = exp( gam1 + gam2*x(i) );
        v(i) = randn(1,1)*sqrt(sig2);
        u(i) = rho*u(i-1) + v(i);
        y(i) = beta1 + beta2*x(i) + u(i);

    end

    % Estimate the model by OLS
    X     = [ ones(t,1) x ];
    bols  = X\y;
    e     = y - X*bols;
    s2    = mean(e.^2);
    omols = s2*inv(X'*X);
    seols = sqrt(diag(omols)); 


    % Estimate the model by ML
    bstart = [ bols ; log( s2 ) ];

    % Optimization settings
    options          = optimset( 'LargeScale','off','Display','off');
    [bhat,~,~,~,~,H] = fminunc( @(b) neglog(b,y,x),bstart,options);
    
    % Computing Standard Errors
    % Based on the Hessian
    iH  = inv( H );
    seH = sqrt( (1/t)*diag( iH ) );

    
    % Based on the OPG matrix
    g   = numgrad( @lnlt,bhat,y,x );
    j   = g'*g/t;
    seG = sqrt( (1/t)*diag(inv( j )) );

    % Based on QMLE
    seQ = sqrt( (1/t)*diag( iH*j*iH ) );
    
    disp('    Standard Errors ');
    disp('    Hessian   OPG       QMLE');
    disp('------------------------------');
    disp( [ seH seG seQ ] );

    
    %  Compute qmle standard errors with lag length based on p 
    p = floor( 4*(t/100)^(2/9) );
 
    for i = 1:p;
 
          gmat = g((i+1):t,:)'*g(1:(t-i),:)/t;
          j    = j + (1.0 - i/(p+1))*(gmat + gmat');
    
    end
    seQp = sqrt( (1/t)*diag( iH*j*iH) );
     
    disp('Standard errors (QMLE - p)');
    disp('-------------------------');
    disp( seQp )

    disp('Standard errors OLS');
    disp('-------------------------');
    disp( seols )

    
end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
% Negative unconstrained log-likelihood
%-------------------------------------------------------------------------
function logl = neglog( b,y,x )

    logl = -mean( lnlt( b,y,x ) ); 
    
end

%-------------------------------------------------------------------------
%   Log-likelihood at every observation
%-------------------------------------------------------------------------
function loglt = lnlt( b,y,x )
    
     u     = y - b(1) - b(2)*x;
     sig2  = exp( b(3) ); 
     z     = u./sqrt( sig2 );
     loglt = - 0.5*log( 2*pi ) - 0.5*log( sig2 ) - 0.5*z.^2 ;
    
end



