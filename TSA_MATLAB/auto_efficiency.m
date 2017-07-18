%=========================================================================
%
%  Simulation example to compare the asymptotic distribution of the MLE 
%  and the OLS estimators of the autocorrelated regression model for 
%  alternative assumptions about the explanatory variable xt
%
%=========================================================================
function auto_efficiency( )

    clear all;
    clc;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',1) )

    %  Simulate a regression model with an AR(1) disturbance term      
    t = 200;                                                    
    beta0  = 1.0;
    beta1  = 1.0;
    rho1   = 0.6;
    sigma2 = 0.0036;
    theta0 = [beta0; beta1; rho1; sigma2];        
    phi1   = 0.6;
    sigma2_w = 0.0036;

    % Two instances where xt is treated as fixed
    %x = (1:1:t)'%.*randn(t,1);  
    % x = sin(2*pi*[0:1:t-1]'/t);     

    ndraws = 500;              
    theta_exact  = zeros(ndraws,4);
    theta_cond   = zeros(ndraws,4);
    theta_ols    = zeros(ndraws,4);

    for k = 1:ndraws
    
        %   Generate data       
        v = sqrt(sigma2)*randn(t,1);           
        u = recserar( v , sqrt(1/(1-rho1^2))*v(1) , rho1 );
        w  = sqrt(sigma2_w)*randn(t,1);
        
        % For the case where xt is not fixed
        x = recserar( w , sqrt(1/(1-phi1^2))*w(1) , phi1 );      

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

    %  Compute statistics of sampling distributions and print results    
    mse_exact = mean(((theta_exact - repmat(theta0',ndraws,1)).^2));
    mse_cond  = mean(((theta_cond  - repmat(theta0',ndraws,1)).^2));
    mse_ols   = mean(((theta_ols   - repmat(theta0',ndraws,1)).^2));

disp('Beta0       Beta1        Rho1       Sigma^2');
disp(['Population parameter       = ' num2str(theta0')]);
disp('--------------------------------------------------------------------');
disp(['Mean    (exact MLE)        = ' num2str(mean(theta_exact))]);
disp(['Bias (exact MLE)           = ' num2str(100*(mean(theta_exact)-theta0'))]);
disp(['MSE  (exact MLE)           = ' num2str(100*mse_exact)]);
disp(['RMSE    (exact MLE)        = ' num2str(sqrt(mse_exact))]);
disp(' ');
disp(['Mean    (cond. MLE)        = ' num2str(mean(theta_cond))]);
disp(['Bias (cond. MLE)           = ' num2str(100*(mean(theta_cond)-theta0'))]);
disp(['MSE  (cond. MLE)           = ' num2str(100*mse_cond)]);
disp(['RMSE    (cond. MLE)        = ' num2str(sqrt(mse_cond))]);
disp(' ');
disp(['Mean (OLS)                 = ' num2str(mean(theta_ols))]);
disp(['Bias (OLS)                 = ' num2str(100*(mean(theta_ols)-theta0'))]);
disp(['MSE  (OLS)                 = ' num2str(100*mse_ols)]);
disp(['RMSE (OLS)                 = ' num2str(sqrt(mse_ols))]);
disp(' ');
disp(['Efficiency (cond/exact)    = ' num2str((mse_cond./mse_exact))]);
disp(['Efficiency (ols/exact)     = ' num2str((mse_ols./mse_exact))]);

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