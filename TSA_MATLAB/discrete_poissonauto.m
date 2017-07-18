%=========================================================================
%
%   Program to generate the finite sampling distribution of the 
%   maximum likelihood estimator of the Poisson autoregressive model 
%   using the binomial thinning operator.
%
%   Results also given for the conditional least squares estimator.
%
%=========================================================================
function discrete_poissonauto( )

    clear all
    clc
    
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',12) );

    % Parameters
    t      = 100;                                                                      
    ndraws = 500;                          
    rho    = 0.3;                              
    lam    = 3.5;
    theta0 = [rho ; lam ];                          

    % Loop to generate sampling distributions   
    theta_mle = zeros(ndraws,2);
    theta_cls = zeros(ndraws,2);

    % Starting values (transformed to satisfy domain restrictions)
    theta_0 = [norminv(rho,0,1) ; log(lam) ];      

    h = waitbar(0,'Please wait...');
    
    for j = 1:ndraws

        waitbar(j/ndraws)
        % Generate data
        u = poissrnd(lam,t,1);     

        % Initialize y by choosing the median of draws from the 
        % unconditional distribution to ensure that the starting value is positive
        y = ones(t,1)*median(poissrnd(lam/(1-rho),t,1));     

        y = uint8( y );
        
        i = 2;
        while i <= t

            e = rand(y(i-1),1);   
            b_thin = sum( e < rho );                          
            y(i) = b_thin + u(i);      

            % Ensure that y(i) is a positive 
            if y(i) > 0
                i = i + 1;               
            end
            
        end
        
        % Estimate by MLE
        ops = optimset('LargeScale','off','Display','off');
    
        theta          = fminsearch(@(b) neglog(b,y),theta_0,ops);
        theta_mle(j,:) = [ normcdf(theta(1),0,1)  exp(theta(2)) ];                       


        % Estimate by CLS     
        xx = [ double(trimr(y,0,1))   ones(t-1,1) ] ;
        yy = double(trimr(y,1,0));
      
        b_cls = xx\yy;     
        theta_cls(j,:) = b_cls';
        
    end
    close( h )
    
    % Compute statistics of sampling distributions  
    mse_mle    = mean( bsxfun(@minus,theta_mle,theta0').^2 );
    rmse_mle   = sqrt( mse_mle );
    mse_cls    = mean( bsxfun(@minus,theta_cls,theta0').^2 );
    rmse_cls   = sqrt( mse_cls );

    disp(' ' );
    disp(['MSE  (MLE)                 = ',num2str(mse_mle) ]);
    disp(['RMSE (MLE)                 = ',num2str(rmse_mle) ])

    disp(' ' );
    disp(['MSE  (CLS)                 = ',num2str(mse_cls) ]);
    disp(['RMSE (CLS)                 = ',num2str(rmse_cls) ])

    disp(' ' );
    disp(['Efficiency (mle/cls)    = ',num2str((mse_mle./mse_cls))]);

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Log-likelihood function of Poisson autoregression
%-------------------------------------------------------------------------
function lf = neglog(b,y)

    y = double(y);
    % Restrict domain of parametes
    rho = normcdf(b(1),0,1);   
    lam = exp(b(2));           

    t = length(y);
    f = zeros(t,1);
 
    for i = 2:t           

        sum = 0.0;      
        tmp = [ y(i); y(i-1) ];
        
        for k = 0:min(tmp)
            
            sum1 = factorial(y(i-1))/(factorial(k)*factorial(y(i-1)-k));
            sum2 = rho^k*(1-rho)^(y(i-1)-k);
            sum3 = lam^(y(i)-k)*exp(-lam)/factorial(y(i)-k);
            sum  = sum + sum1*sum2*sum3;
        end

        f(i) = sum;

    end

    % Exclude first observation from the likelihood
    f  = trimr(f,1,0);               
   	lf = -mean( log(f) );

end

