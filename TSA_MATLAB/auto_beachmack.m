%=========================================================================
%
%     Simulation example to reproduce the Beach and MacKinnon (1978)
%     Econometrica, pp.51-58 study which derives the sampling
%     distributions MLE estimators of regression models with
%     autocorrelation.
%
%=========================================================================

function auto_beachmack( )

    clear all;
    clc;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',12345) )

    % Simulate a regression model with an AR(1) disturbance term     
    t = 20;                                                                                        
    beta0  = 1.0;
    beta1  = 1.0;
    rho1   = 0.6;
    sigma2 = 0.0036;
    theta0= [beta0 ; beta1 ; rho1 ; sigma2];      
    
    % xt is fixed in repeated samples
    x = exp(0.04*(1:1:t)') + sqrt(0.0009)*randn(t,1);   
                                                     

    ndraws = 200;          % used by Beach and MacKinnon     
    theta_exact = zeros(ndraws,4);
    theta_cond  = zeros(ndraws,4);
    theta_ols   = zeros(ndraws,4);

    for k = 1:ndraws

        %     Generate data      
        v = sqrt(sigma2)*randn(t,1);                                                                                                             
        u = recserar( v , sqrt(1/(1-rho1^2))*v(1) , rho1 );                
        y = beta0 + beta1*x + u;

        % Exact MLE
        flag = 0;        
        theta = fminsearch(@(b) neglog(b,y,x,flag),theta0);       
        theta_exact(k,:) = [ theta(1:2) ; tanh(theta(3)) ; abs(theta(4))]';       

        % Conditional MLE       
        flag =  1;
        theta = fminsearch(@(b) neglog(b,y,x,flag),theta0);       
        theta_cond(k,:) = [ theta(1:2) ; tanh(theta(3)) ; abs(theta(4)) ]'; 
        
        % OLS       
        b_ols = [ones(t,1),x]\y;
        e    = y - [ones(t,1),x]*b_ols;
        sig2_ols = mean(e.^2);
        rho_ols = e(1:end-1)\e(2:end);
        theta_ols(k,:) = [b_ols' , sig2_ols , rho_ols];
      
    end

    % Compute statistics of sampling distributions and print results     
    rmse_exact = sqrt( mean( (mean(theta_exact) - theta0').^2 ) );
    rmse_cond  = sqrt( mean( (mean(theta_cond) - theta0').^2 ) );
    rmse_ols   = sqrt( mean( (mean(theta_ols) - theta0').^2 ) );


    fprintf ('                             Beta0      Beta2      Rho1      Sigma^2\n');
    fprintf ('Population parameter        %2.3f      %2.3f      %2.3f     %2.3f\n\n\n', theta0');
    fprintf ('Mean (exact MLE)        %2.3f      %2.3f      %2.3f     %2.3f\n', mean(theta_exact)');
    fprintf ('Bias (exact MLE)         %2.3f      %2.3f      %2.3f     %2.3f\n', theta0'-mean(theta_exact));
    disp(' ')
    fprintf ('Mean (cond. MLE)         %2.3f      %2.3f      %2.3f     %2.3f \n', mean(theta_cond)');
    fprintf ('Bias (cond. MLE)         %2.3f      %2.3f      %2.3f     %2.3f \n', theta0'-mean(theta_cond));
    disp(' ');
    fprintf ('Mean (OLS)               %2.3f      %2.3f      %2.3f     %2.3f \n', mean(theta_ols)');
    fprintf ('Bias (OLS)               %2.3f      %2.3f      %2.3f     %2.3f \n', theta0'-mean(theta_ols));
    disp(' ');
    fprintf ('\nRMSE (exact MLE)        %f \n', rmse_exact');
    fprintf ('\nRMSE (cond. MLE)        %f \n', rmse_cond');
    fprintf ('\nRMSE (OLS)              %f \n\n', rmse_ols');
    disp(' ')
    fprintf ('Efficiency (cond/exact) %f \n', (rmse_cond./rmse_exact)');
    fprintf ('Efficiency (ols/exact)  %f \n', (rmse_ols./rmse_exact)');

end
%
%-------------------------Functions------------------------------------
%
%-----------------------------------------------------------------------
%      Negative log-likelihood function   
%-----------------------------------------------------------------------

function lnl = neglog(b,y,x,flag)

     beta0 = b(1);                       
     beta1 = b(2);
     rho1  = tanh(b(3));   %  rho1 stays in the unit circle
     sig2  = abs(b(4));    %  variance is positive     

     u     = y - beta0 - beta1*x;
     v     = u(2:end) - rho1*u(1:end-1);
     lnl_0 = - 0.5*log(2*pi) - 0.5*log(sig2) + 0.5*log(1 - rho1^2) ...
             - 0.5*(u(1) - 0).^2/(sig2/(1 - rho1^2));                
     lnl_1 = - 0.5*log(2*pi) - 0.5*log(sig2) - 0.5*v.^2/sig2;                                                           

     if flag
        lnl = -mean(lnl_1);
     else
        lnl = -mean( [lnl_0 ; lnl_1] );
     end

end