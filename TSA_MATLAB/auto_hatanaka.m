%=========================================================================
%  
%  Simulation example to compure the asymptotic distribution of the MLE
%  estimator and the Hatanaka estimator for the regression model with 
%  autocorrelation.
%
%=========================================================================

function auto_hatanaka( )

    clear all;
    clc;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',1) )

    t = 1000000;     

    beta0  = 1.0;
    beta1  = 1.0;
    alpha1 = 0.5;
    rho1   = 0.6;
    sigma2 = 0.1;
    theta0 = [beta0, beta1, alpha1, rho1, sigma2];                                                            
    
    x = rand(t,1) - 0.5;                
    v = sqrt(sigma2)*randn(t,1);    
    u = recserar(v, sqrt(1/(1-rho1^2))*v(1), rho1);  % Beach and MacKinnon, fn3     
    y = recserar(beta0 + beta1*x + u , 0.0, alpha1);                                                                      
    cov_theoretical = rho1*sigma2/((1-alpha1*rho1)*(1-rho1^2));
    cov_simulated   = cov(trimr(y,0,1), trimr(u,1,0)); 

    disp('Theoretical and simulated covariance');
    disp([ cov_theoretical, cov_simulated(1,2) ]);

    %   Monte Carlo replications  
    t      = 1000;
    ndraws = 1000;  
    x      = rand(t,1) - 0.5;
    theta_cond     = zeros(ndraws,5);
    theta_hatanaka = zeros(ndraws,5);
    theta_ols      = zeros(ndraws,5);

    for k = 1:ndraws

        %   Generate data       
        v = sqrt(sigma2)*randn(t,1);                                                                                                          
        u = recserar( v , sqrt(1/(1-rho1^2))*v(1) , rho1 );   
        y = recserar( beta0 + beta1*x + u , 0.0, alpha1 );    
        
        %   Conditional MLE      
        theta = fminsearch(@(b) neglog(b,y,x),theta0);
        theta_cond(k,:) = [theta(1) theta(2) tanh(theta(3)) tanh(theta(4))  abs(theta(5))];                         

        %   Hatanaka 2-step efficient estimator       
        theta = hatanaka(y,x,t);
        theta_hatanaka(k,:) = theta;

        %   OLS       
        b_ols = [ones(t-1,1), trimr(x,1,0), trimr(y,0,1)]\trimr(y,1,0);
        u_hat = trimr(y,1,0) - [ones(t-1,1), trimr(x,1,0), trimr(y,0,1)]*b_ols;
        rho_ols = trimr(u_hat,0,1)\trimr(u_hat,1,0);
        sig2_ols = mean(u_hat.^2);
        theta_ols(k,:) = [b_ols', rho_ols, sig2_ols];

    end

    % Compute statistics of sampling distributions and print results     
    mse_cond     = mean((theta_cond - repmat(theta0,ndraws,1)).^2);
    mse_hatanaka = mean((theta_hatanaka - repmat(theta0,ndraws,1)).^2);
    mse_ols      = mean((theta_ols - repmat(theta0,ndraws,1)).^2);

    disp('              Beta0      Beta1      Alpha1     Rho1      Sigma^2');
    disp(['Population parameter            ' num2str(theta0)]);
    disp(['Mean       (cond. MLE)          ' num2str(mean(theta_cond))]);
    disp(['Bias(x100) (cond. MLE)          ' num2str(100*(mean(theta_cond)-theta0))]);
    disp(['MSE(x100)  (cond. MLE)          ' num2str(100*mse_cond)]);
    disp(['RMSE(x100) (cond. MLE)          ' num2str(100*sqrt(mse_cond))]);
    
    disp(['Mean       (Hatanaka)           ' num2str(mean(theta_hatanaka))]);
    disp(['Bias(x100) (Hatanaka)           ' num2str(100*(mean(theta_hatanaka)-theta0))]);
    disp(['MSE(x100)  (Hatanaka)           ' num2str(100*mse_hatanaka)]);
    disp(['RMSE(x100) (Hatanaka)           ' num2str(100*sqrt(mse_hatanaka))]);
    disp(['Mean (OLS)                      ' num2str(mean(theta_ols))]);
    disp(['Bias(x100) (OLS)                ' num2str(100*(mean(theta_ols)-theta0))]);
    disp(['MSE(x100)  (OLS)                ' num2str(100*mse_ols)]);
    disp(['RMSE(x100) (OLS)                ' num2str(100*sqrt(mse_ols))]);
    disp(['Efficiency (hatanaka/cond)      ' num2str((sqrt(mse_hatanaka)./sqrt(mse_cond)))]);
    disp(['Efficiency (ols/cond)           ' num2str((sqrt(mse_ols)./sqrt(mse_cond)))]);

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%   Log-likelihood function at each observation            
%-------------------------------------------------------------------------
function lf = neglog(b,y,x)
        
    beta0  = b(1);                        
    beta1  = b(2);
    alpha1 = tanh(b(3));     % rho1 stays in the unit circle    
    rho1   = tanh(b(4));     % rho2 stays in the unit circle    
    sig2   = abs(b(5));      % variance is positive     
    u      = trimr(y,1,0) - beta0 - beta1*trimr(x,1,0) - alpha1*trimr(y,0,1);
    v      = trimr(u,1,0) - rho1*trimr(u,0,1);
    lnl    = - 0.5*log(2*pi) - 0.5*log(sig2) - 0.5*v.^2/sig2;
    lf     = -mean( lnl );
     
end
%-------------------------------------------------------------------------
% Hatanaka 2-step efficient estimator      
%-------------------------------------------------------------------------
function h1 = hatanaka(y,x,t)
    
    yvar = trimr(y,1,0);
    xvar = [ones(t-1,1), trimr(x,1,0), trimr(y,0,1)];
    zvar = [ones(t-1,1), trimr(x,1,0), trimr(x,0,1)];
    
    %   IV initial estimates of mean parameters
    biv  = inv(zvar'*xvar)*(zvar'*yvar);         
    u    = yvar - xvar*biv;
    
    % Estimate rho
    rho  = trimr(u,0,1)\trimr(u,1,0);     
    yvar = trimr(y,2,0) - rho.*trimr(y,1,1);             
    xvar = [ones(t-2,1), (trimr(x,2,0) - rho.*trimr(x,1,1)), (trimr(y,1,1) - rho.*trimr(y,0,2)), trimr(u,0,1)];
    
    % Regression on transformed variables and update rho
    b    = xvar\yvar;                         
    rho  = rho + b(4);                           

    %    Compute residual variance 
    v    = trimr(y,2,0)-([ones(t-2,1), trimr(x,2,0), trimr(y,1,1)]*b(1:3) -rho*(trimr(y,1,1)-[ones(t-2,1), trimr(x,1,1), trimr(y,0,2)]*b(1:3)));
    sig2 = mean(v.^2);                                    
    h1   = [(b(1:3))', rho, sig2];

end