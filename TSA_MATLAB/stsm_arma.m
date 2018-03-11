%=========================================================================
%
%   Program to estimate an ARMA model
%
%=========================================================================
function stsm_arma( )

    clear all
    clc

    % Read the data: quarterly US data from Jan-1959 to Dec-1998
    load sims_data.mat
 
    % Define variables
    r    = ytdata(:,1);        
    lex  = log( ytdata(:,2) );
    lcp  = log( ytdata(:,3) );
    lm   = log( ytdata(:,4) );
    lp   = log( ytdata(:,5) );
    lo   = log( ytdata(:,6) );
    sdum = ytdata(:,7:17);
    
    t = length(ytdata);

    % Choose y variable
    y = r;                                          % Interest rate
    %y = 100*(trimr(lo,12,0) - trimr(lo,0,12));      % Annual output growth
    %y = 100*(trimr(lp,12,0) - trimr(lp,0,12));      % Annual inflation 
        
    
    bstart = 0.1*ones(3,1);
    options = optimset('LargeScale','off','Display','iter');
    
    [ theta,fval,~,~,~,H ]  = fminunc( @(b) neglogl(b,y),bstart,options );  % H  «Hessianæÿ’Û 
    
    lnl1 = -fval;
    vc   = (1/(t-1))*inv(H);

    % Wald test (AR)
    r = [ 0  1  0 ];
    q = 0;
    w = (r*theta - q)'*inv(r*vc*r')*(r*theta - q);

    disp(' ')
    disp(['Wald statistic (AR)     = ', num2str(w) ]);
    disp(['p-value                 = ', num2str(1-chi2cdf(w,1)) ]);
    disp(' ')

    % Wald test (MA)   
    r = [ 0  0  1 ];
    q = 0;
    w = (r*theta - q)'*inv(r*vc*r')*(r*theta - q);
    disp(' ')
    disp(['Wald statistic (MA)     = ', num2str(w) ]);
    disp(['p-value                 = ', num2str(1-chi2cdf(w,1)) ]);
    disp(' ')
 
    % Wald test (Joint)  
    r = [ 0  1  0 ;
          0  0  1 ];
    q = [ 0 ; 0 ];
    w = (r*theta - q)'*inv(r*vc*r')*(r*theta - q);
    disp(' ')
    disp(['Wald statistic (Joint)  = ', num2str(w) ]);
    disp(['p-value                 = ', num2str(1-chi2cdf(w,2)) ]);
    disp(' ')

    % Likelihood Ratio test  
    lnl0 = -neglogl( [mean(y) 0 0],y ); 
    lr   = -2*(t-1)*(lnl0 - lnl1);

    disp(' ')
    disp(['LR statistic (Joint)    = ', num2str(lr) ]);
    disp(['p-value                 = ', num2str(1-chi2cdf(lr,2)) ]);
    disp(' ')
    
end
%--------------------------- Subroutines ----------------------------------
% 
%--------------------------------------------------------------------------
% Log-likelihood function for an ARMA model
%--------------------------------------------------------------------------

function f = neglogl( b,y )

    [ t,n ] = size( y );
    v       = zeros( t,1 );
    lf      = zeros( t-1,1 );
    
    % First loop over MA term   
    for i = 2:t
        
        v(i) = y(i)-b(1)-b(2)*y(i-1)-b(3)*v(i-1);
    end
    
    v       = trimr( v,1,0 );
    vc      = v'*v/(t-1);

    for i = 1:t-1;

        lf(i) = -0.5*n*log(2*pi) - 0.5*log(det(vc)) ...
                - 0.5*v(i,:)*inv(vc)*v(i,:)';    
    end
    
    f = -mean( lf );
end
