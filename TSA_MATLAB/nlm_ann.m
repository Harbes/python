%=========================================================================
%
%      Estimate an artificial neural network by maximum likelihood  
%
%=========================================================================
function nlm_ann( )

    clear all;
    clc;

    % Initialise the random number generator  
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234) ); 	
    
    % Parameters
    nobs = 2000;                         
    phi    = 0.2;                    
    gam    = 2.0;
    delta0 = -2.0;
    delta1 = 2.0;
    sig2   = 0.1;


    % Simulate  data						     
    y = zeros(nobs+100,1) ;
    w = zeros(nobs+100,1) ;
    u = sqrt(sig2)*randn(nobs+100,1);

    for t = 2:nobs+100

        w(t-1) = 1/( 1 + exp( - ( delta0 + delta1*y(t-1) ) ) );
        y(t)   = phi*y(t-1) + gam*w(t-1) + u(t);
    end
    y = trimr(y,100,0);

    % Estimate the model using BFGS
    theta_0 = [phi ; gam ; delta0 ; delta1];
    ops     = optimset( 'LargeScale','off','Display','off');  
    
    [theta,lnl] = fminunc(@(b) neglog(b,y),theta_0,ops);
    
    disp('    True    Estimated' )
    disp([theta_0 theta]);
    disp(-lnl);

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Log-likelihood function
%-------------------------------------------------------------------------
function lf = neglog(b,y)

    nobs = length(y);
    v    = zeros(nobs,1);
    w    = zeros(nobs,1);

    for t = 2:nobs
         w(t-1) = 1/( 1 + exp( - ( b(3) + b(4)*y(t-1) ) ) );
         v(t)   = y(t) - b(1)*y(t-1) - b(2)*w(t-1);
    end
     v  = trimr(v,1,0);	 
     s2 = mean(v.^2)'; 
     z  = v./sqrt(s2);
     lf = -mean( -0.5*log(2*pi)-0.5*log(s2)-0.5*z.^2 );

end
     


