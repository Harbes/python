%=========================================================================
%
%   Recursions of the univariate Kalman filter with smoothing
%
%=========================================================================
function lfac_smooth( )

    clear all
    clc

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234567) );

    % Simulate the data       
    t    = 5;
    lam  = 1.0;
    sig  = 0.5;
    phi  = 0.8;
    sigv = 1;

    u = sig*randn(t,1);
    v = randn(t,1);
    f = recserar(v,0.0,phi);
    y = lam*f + u;

    % Reproduce the GAUSS results in text
    %y =  [   0.079 -0.091 -0.407 2.070 -0.047 ];

    
    % Return the smoothed latent factor     
    theta0 = [lam ; sig ; phi ; sigv];
    lfac = lnlt(theta0,y);
    
    figure;
    plot([y lfac])

end
%--------------------------- Functions ----------------------------------
% 
%--------------------------------------------------------------------------
% Univariate Kalman filter
%--------------------------------------------------------------------------
function lfac = lnlt(b,y)

    % Unpack parameters
	lam  = b(1);
    sig  = b(2)^2;
    phi  = b(3);
    sigv = b(4)^2;    
    
    %str = ['Actual   ';'Predict  ';'Update   ';'Smooth   '];
        
    % Allocate arrays
    t   = length(y);
    lnl = zeros(t,1);
    spr = zeros(t,1);     %     Predicted factor estimates     
    sup = zeros(t,1);     %     Updated factor estimates       
    ssm = zeros(t,1);     %     Smoothed factor estimates      
    ppr = zeros(t,1);     %     Predicted factor cov  
    pup = zeros(t,1);     %     Updated factor cov 
    psm= zeros(t,1);      %     Smoothed factor cov 


    % Recursions of the Kalman Filter
    % Initialisation following Harvey and Hamilton
    st     = 0.0;
    pt     = 1/(1-phi^2);   
    spr(1) = st;
    ppr(1) = pt;   
    mt     = lam*st;
    vt     = lam*pt*lam + sig;
    ut     = y(1) - mt;
    lnl(1) = - 0.5*log(2*pi) - 0.5*log(det(vt)) - 0.5*ut^2/vt;
    kal    = pt*(lam/vt);
    s0     = st + kal*ut;
    p0     = pt - kal*lam*pt;
    sup(1) = s0;
    pup(1) = p0;
       
    
    % Main loop over observations
    for i = 2:t
    
        % Predict and store predictions
        st      = phi*s0;
        pt     = phi*p0*phi + sigv;
        spr(i) = st;
        ppr(i) = pt;
                  
        % Observation     
        mt = lam*st;
        vt = lam*pt*lam + sig;
        ut = y(i) - mt;
        
        % Log-likelihood function
        lnl(i) = - 0.5*log(2*pi) - 0.5*log(det(vt)) - 0.5*ut^2/vt;

                  
        % Update and store updates
        kal = pt*lam/vt;
        s0 = st + kal*ut;
        p0 = pt - kal*lam*pt;
        sup(i) = s0;
        pup(i) = p0;

    end
    
    
    % Backward recursion to smooth and extract the factor
    ssm(t) = sup(t);   
    for i = t-1:-1:1        
    
        j      = (pup(i)*phi)/ppr(i+1);       
        ssm(i) = sup(i) + j*(ssm(i+1) - spr(i+1))*j;
        
    end
    
    lfac = ssm;

end

