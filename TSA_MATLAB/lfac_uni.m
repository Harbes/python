%=========================================================================
%
%   Recursions of the univariate Kalman filter
%
%=========================================================================
function lfac_uni( )

    clear all
    clc

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234567) );

    % Simulate the data       
    t    = 5;
    lam  = 1.0;
    sig  = 0.5;
    phi  = 0.8;
    sige = 1;

    u = sig*randn(t,1);
    v = randn(t,1);
    f = recserar(v,0.0,phi);
    y = lam*f + u;

    % Reproduce the GAUSS results  
    %y =  [   0.079 -0.091 -0.407 2.070 -0.047 ];

    
    % Kalman filter     
    theta0 = [lam ; sig ; phi ; sige];
    lnlt(theta0,y);

end
%--------------------------- Functions ----------------------------------
% 
%--------------------------------------------------------------------------
% Univariate Kalman filter
%--------------------------------------------------------------------------
function lnlt(b,y)

    % Unpack parameters
	lam  = b(1);
    sig  = b(2)^2;
    phi  = b(3);
    sigv = b(4)^2;    
    
    str = ['yt    ';'st    ';'pt    ';'mt    ';'vt    ';'ut    ';'kal   ';'lnlt  '];
        
    % Allocate arrays
    t   = length(y);
    lnl = zeros(t,1);

    % Recursions of the Kalman Filter
    % Initialisation following Harvey and Hamilton
    st     = 0.0;
    pt     = 1/(1-phi^2);   
    mt     = lam*st;
    vt     = lam*pt*lam + sig;
    ut     = y(1) - mt;
    lnl(1) = - 0.5*log(2*pi) - 0.5*log(det(vt)) - 0.5*ut^2/vt;
    kal    = pt*(lam/vt);
    s0     = st + kal*ut;
    p0     = pt - kal*lam*pt;
       
    % Display initialisation
    disp(' ');
    disp('Filter initialisation');
    tmp = [y(1) st pt mt vt ut kal lnl(1)]';
    disp([ str num2str(tmp) ]);
    disp(' ');
    
    % Main loop over observations
    for i = 2:t
    
        % Prediction  
        st = phi*s0;
        pt = phi*p0*phi + sigv;
                  
        % Observation     
        mt = lam*st;
        vt = lam*pt*lam + sig;
        ut = y(i) - mt;
        
        % Log-likelihood function
        lnl(i) = - 0.5*log(2*pi) - 0.5*log(det(vt)) - 0.5*ut^2/vt;

                  
        % Updating
        kal = pt*lam/vt;
        s0 = st + kal*ut;
        p0 = pt - kal*lam*pt;

        % Display results
        disp(' ');
        disp(['Iteration  = ',num2str(i) ]);
        tmp = [y(i) st pt mt vt ut kal lnl(i)]';
        disp([ str num2str(tmp) ]);
        disp(' ');
    end
end

