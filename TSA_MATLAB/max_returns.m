%=========================================================================
%
%   Program to Program to estimate the distribution of asset returns 
%   using a standardised Student t distribution
%
%   Asset price data from 6 August 2010 to 2 January 2001 (note that the
%   data are in reverse order ie from recent to past)
%
%=========================================================================
function max_returns( )

    clear all
    clc

    % Load data
    load diversify.mat

    t = 2413;
    
    % Select appropriate sample  
    pt_apple     = pt_apple(1:t);
    pt_ford      = pt_ford(1:t);

    % Compute percentage returns  
    r_apple = 100*diff(log(pt_apple)); 
    r_ford  = 100*diff(log(pt_ford));

    % Select a return 
    y = r_apple;

    % Compute desrciptive statistics   
    m = mean(y);
    s = std(y);
    z = (y - m)/s;

    skew = mean( z.^3);
    kurt = mean( 4.^3);

    disp(['Mean     =   ',num2str(m) ]);
    disp(['Std dev. =   ',num2str(s) ]);
    disp(['Skewness =   ',num2str(skew) ]);
    disp(['Kurtosis =   ',num2str(kurt) ]);
    disp(['Maximum  =   ',num2str(max(z)) ]);
    disp(['Minimum  =   ',num2str(max(z)) ]);

    figure(1)    
    histfit( z, 21);


    % Estimate the model by MLE with Student t distribution 
    start = [m ; s ; 5 ];
    ops   = optimset('LargeScale','off','Display','off');
    [ bhat,lf,a,a,a,hess] = fminunc(@(b) lnl( b,z ),start,ops );

    lf = -lf;
    vc = 1/(t-1)*inv(hess);
    
    disp('Unrestricted model results')
    disp(['Log-likelihood function = ',num2str(lf) ]);
    disp('    Param     Std err.')
    disp( [ bhat sqrt(diag(vc))] )
    
    % Estimate the model by MLE with normal distribution 
    [ bhat,lf,a,a,a,hess] = fminunc(@(b) lnl1( b,z ),start(1:2),ops );

    lf = -lf;
    vc = 1/(t-1)*inv(hess);
    
  
    disp('Restricted model results')
    disp(['Log-likelihood function = ',num2str(lf) ]);
    disp('    Param     Std err.')
    disp( [ bhat sqrt(diag(vc))] )

end
% ------------------------ Functions ------------------------------------%
%
%-------------------------------------------------------------------------
%  Log-likelihood of a Student t distribution
%-------------------------------------------------------------------------
function lf = lnl(b,z)

    
        m   = b(1);                  
        s   = abs(b(2));              
        v   = abs(b(3));              
        s2  = s^2;                   
        tmp = log( gamma((v+1)/2) ) - 0.5*log(pi) - log(gamma(v/2)) - 0.5*log(v-2) ...
            - 0.5*log(s2) - ((v+1)/2)*log( 1 + ((z - m).^2)./(s2*(v-2)) );

        lf  = -mean(tmp);

end

%-------------------------------------------------------------------------
%  Log-likelihood of a normal distribution
%-------------------------------------------------------------------------
function lf = lnl1(b,z)

    
        m   = b(1);                  
        s   = abs(b(2));                                           
        tmp = log(normpdf(z,m,s));

        lf  = -mean(tmp);

end
