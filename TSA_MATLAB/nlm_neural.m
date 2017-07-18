%=========================================================================
%
%      Estimate an artificial neural network using the neural algorithm  
%
%=========================================================================
function nlm_neural( )

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

    % Estimate the model using neural algorithm
    n    = 10000;                                     
    pmat = zeros(n,5);          
    ylag = trimr(y,0,1);         
    y    = trimr(y,1,0);         


    for iter = 1:n
        
        d0 = -5 + 10*rand(1,1);
        d1 = -5 + 10*rand(1,1);
        w  = 1./( 1 + exp( - ( d0 + d1*ylag ) ) );   
        x  = [ylag   w ];
        b  = x\y;
        v  = y - x*b;
        s2 = v'*v/length(y);

        pmat(iter,:) = [b'   d0   d1   s2 ];
    end

    % Sort the parameter estimates based on s2    
    pmat = sortrows(pmat,5); 

    disp('Parameter estimates')
    disp(['phi         = ',num2str(pmat(1,1))]);
    disp(['gam         = ',num2str(pmat(1,2))]);
    disp(['delta0      = ',num2str(pmat(1,3))]);
    disp(['delta1      = ',num2str(pmat(1,4))]);
    disp(['sigma^2     = ',num2str(pmat(1,5))]);
end