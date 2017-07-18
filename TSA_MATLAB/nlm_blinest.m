%=========================================================================
%
%   Estimate and test a bilinear time series model by maximum likelihood
%
%=========================================================================

function nlm_blinest( )

    clear all
    clc
    
    % Initialise the random number generator  
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123457) ); 	

    % Simulate data
    t = 1000;
    y = zeros(t+100,1) ;
    v = sqrt(0.1)*randn(t+100,1);
    
    for k = 2:t+100
 
        y(k) = 0.1 + 0.4*y(k-1) + 0.4*y(k-1)*v(k-1) + v(k);
    end
    y = trimr(y,100,0);

    %  Estimate the unrestricted model 
    ops     = optimset('LargeScale','off','Display','off');
    theta_0 = [0.1 ; 0.1 ; 0.1];
    [theta1,lf1,~,~,~,hess] = fminunc(@(b) neglog(b,y),theta_0,ops); 

    vc  = (1/t)*inv(hess);
    se  = sqrt(diag(vc));
    lf1 = -lf1;
    
    disp(['Unrestricted log-likelihood     = ',num2str(lf1) ])
    disp('  Estimates     Std.Errors ' )
    disp([ theta1 se ])
    
    % Restricted model
    theta_0 = [0.1 ; 0.1 ; 0];

    [~,lf0] = fminunc(@(b) neglogr(b,y),theta_0); 

    lf0 = -lf0;

    % LR test       
    lr = -2*t*(lf0 - lf1);

    disp(['Likelihood ratio test  =  ', num2str(lr) ]);
    disp(['p-value                =  ', num2str(1 - chi2cdf(lr,1)) ]);
    disp(' ');

    % Wald test       
    w = (theta1(3) - 0)^2/vc(3,3);

    disp(['Wald test              =  ', num2str(w) ]);
    disp(['p-value                =  ', num2str(1 - chi2cdf(w,1)) ]);
    disp(' ');

    % LM test       
    % First stage regression      
    x = [ones(size(y,1)-1,1)  trimr(y,0,1)];                     
    y = trimr(y,1,0);           
    k = size(x,2);
    v = y - x*(x\y);                                        

    % Second stage regression      
    x = [trimr(x,1,0)  trimr(x(:,2),1,0).*trimr(v,0,1) ];              
    v = trimr(v,1,0);
    e = v - x*(x\v);                                               

    t  = size(y,1);
    r2 = 1 - sum(e.^2)/sum( (v-mean(v)).^2 );
    lm = t*r2;

    disp(['LM test              =  ', num2str(lm) ]);
    disp(['p-value                =  ', num2str(1 - chi2cdf(lm,1)) ]);

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Unrestricted log-likelihood function
%-------------------------------------------------------------------------
function lf = neglog(b,y)

    nobs = length(y);     
    v    = zeros(nobs,1);
    
    for t = 2:nobs
         v(t) = y(t) - b(1) - b(2)*y(t-1) - b(3)*y(t-1)*v(t-1);
    end
    v = trimr(v,1,0);	 	
    s2 = mean(v.^2)'; 
    z  = v./sqrt(s2);

    lf = -mean( -0.5*log(2*pi)-0.5*log(s2)-0.5*z.^2 );

end
%-------------------------------------------------------------------------
%  Restricted log-likelihood function
%-------------------------------------------------------------------------
function lf = neglogr(b,y)

    nobs = length(y);     
    v    = zeros(nobs,1);
    
    for t = 2:nobs
         v(t) = y(t) - b(1) - b(2)*y(t-1);
    end
    v = trimr(v,1,0);	 	
    s2 = mean(v.^2)'; 
    z  = v./sqrt(s2);

    lf = -mean( -0.5*log(2*pi)-0.5*log(s2)-0.5*z.^2 );

end