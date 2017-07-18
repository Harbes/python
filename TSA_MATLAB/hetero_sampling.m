%=========================================================================
%
%    Simulation example to generate the sampling distributions the MLE 
%    and OLS estimator of a regression model with heteroskedasticity.
%
%=========================================================================
function hetero_sampling( )

    clear all;
    clc;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123456) )

    % Simulate a regression model with multiplicative heteroskedasticity
    t      = 500;                                                                                      
    beta0  = 1.0;
    beta1  = 2.0;
    gam0   = 1.0;
    gam1   = 5.0;
    theta0 = [beta0 ; beta1 ; gam0 ; gam1];            

    % Fixed exogenous variable
    x = rand(t,1);       
                     
    % Monte Carlo replications
    ndraws   = 10000;                
    theta_mle = zeros(ndraws,4);
    theta_ols = zeros(ndraws,4);
    
    % Optimisation settings
    ops = optimset( 'MaxIter',     20000,    ...
                    'MaxFunEvals', 40000);
   
    
    for k = 1:ndraws
        
        % Generate data       
        u = sqrt( exp(gam0 + gam1*x) ).*randn(t,1);                         
        y = beta0 + beta1*x + u;                     	                                   
        
        % ML      
        theta_hat      = fminsearch(@(b) neglog(b,y,x),theta0,ops);       
        theta_mle(k,:) = theta_hat';       
        
        % OLS       
        b_ols          = [ones(t,1),x]\y;
        e              = y - [ones(t,1),x]*b_ols;
        gam0_ols       = log(mean(e.^2));
        theta_ols(k,:) = [b_ols' , gam0_ols , 0];
    end
    
    % Compute statistics of sampling distributions and print results
    rmse_mle   = sqrt( mean( (theta_mle - repmat(theta0',ndraws,1)).^2 ) );
    rmse_ols   = sqrt( mean( (theta_ols - repmat(theta0',ndraws,1)).^2 ) );

    disp(' ML results:')
    disp('    True      Estimate  Bias     RMSE' );
    disp( [ theta0 mean(theta_mle)' mean(theta_mle)'-theta0 rmse_mle' ] );

    disp(' ');
    disp(' OLS results:')
    disp('    True      Estimate  Bias     RMSE' );
    disp( [ theta0 mean(theta_ols)' mean(theta_ols)'-theta0 rmse_ols' ] );
    
    disp(' ');

    disp('Efficiency (OLS/ML)');
    disp((rmse_ols./rmse_mle)');



end
%
%-------------------------Functions------------------------------------
%
%-----------------------------------------------------------------------
%      Negative log-likelihood function   
%-----------------------------------------------------------------------

function lf = neglog( p,y,x )

    lf = -mean( lnlt(p,y,x) );
    
end
%-----------------------------------------------------------------------
%      Log-likelihood function at each observation
%-----------------------------------------------------------------------
function lft = lnlt(p,y,x)

    mu   = p(1) + p(2)*x;
    sig2 = exp(p(3) + p(4)*x);
    lft  = -(0.5)*log(2*pi*sig2) - (y - mu).^2./(2*sig2);        

end



