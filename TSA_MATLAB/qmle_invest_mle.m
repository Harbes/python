%========================================================================
%
%      Example taken from Greene, Table F5.1, 5th edition.
%
%      Estimate the investment model by OLS and compute QMLE standard errors
%
%========================================================================

function qmle_invest_mle( )

% Read in quarterly data for the US over the period 1950II to 2000IV on:

%    Year
%    Quarter
%    Real GDP ($bil)
%    Real consumption expenditures
%    Real investment by private sector
%    Real government expenditures
%    Real disposable personal income
%    Consumer price index
%    Nominal money stock, M1
%    Quarterly average of month end 90 day t bill rate
%    Unemployment rate
%    Population, mil. interpolate of year end figures using constant growth rate per quarter
%    Rate of inflation (first observation is missing)
%    Ex post real interest rate = Tbilrate - Infl. (First observation missing)

%   Table F5.1: Macroeconomics Data Set, Quarterly, 1950I to 2000IV, 204 Quarterly Observations

    clear all
    clc

    load US_MacroData

    ri   = usdata(:,5);
    ry   = usdata(:,3);
    rint = usdata(:,14)/100;
    t    = length(ri);
    trend = (1:1:t)'/100;


    yt = log( ri );
    xt = [ ones( t,1 ) log( ry ) rint trend ];

    
    % Estimate model by OLS
    bols = xt\yt;
    et   = yt - xt*bols;
    s2   = mean(et.^2);

    % Call optimization routine

    bstart = [ -7; 1; -0.5; -0.5; s2 ];

    
    % Optimization settings
    options = optimset( 'LargeScale','off',  ...
                        'Display','final',     ...
                        'MaxIter',15000,     ...
                        'MaxFunEvals',10000, ...
                        'TolFun',1e-12,       ...
                        'TolX',1e-12);
    

%    [bhat,fu] = fminsearch( @(b) lnl(b,yt,xt),bstart,options); 
    [ bhat,fu,flag,out,grad,H ] = fminunc( @(b) lnl( b,yt,xt ),bstart,options);
  
    % Compute standard errors based on the Hessian
    % Use this as an alternative estimate of H
    H = fdhess( @lnl,bhat,yt,xt )


    
    % This is the INVERSE Hessian returned by GAUSS 'maxlik'
%     H = [     1.498473455518  -0.200382212161  -0.044900697510   0.165399548126  -0.000220236028 
%              -0.200382212161   0.026798585773   0.006029688943  -0.022134925958   0.000029452899 
%              -0.044900697510   0.006029688943   0.050609175101  -0.005795199962   0.000006620708 
%               0.165399548126  -0.022134925958  -0.005795199962   0.018402046888  -0.000024322458 
%              -0.000220236028   0.000029452899   0.000006620708  -0.000024322458   0.000000567300  ];
    
    seH = sqrt( diag( inv(t*H) ) );
    
   
    % Compute standard errors based on the OPG matrix
    G   = numgrad( @lnlt,bhat,yt,xt );
    J   = G'*G;
    seG = sqrt( diag( inv( J ) ) );
   
    H = inv( t*H );
    % Compute QMLE standard errors
    
    seQ0 = sqrt( diag( H*J*H) );

 
    %  Compute qmle standard errors with lag length based on p 
    p = floor( 4*(t/100)^(2/9) );
 
    for i = 1:p;
 
          gmat = G((i+1):t,:)'*G(1:(t-i),:);
          J    = J + (1.0 - i/(p+1))*(gmat + gmat');
    
    end
    seQp = sqrt( diag( H*J*H) );
     

    disp('Parameter Estimates');
    disp('-------------------');
    disp( bhat );
 
    
    disp('Standard errors (Hessian)');
    disp('-------------------------');
    disp( seH )
    

    disp('Standard errors (OPG)');
    disp('-------------------------');
    disp( seG )
 
    
    disp('Standard errors (QMLE - p=0)');
    disp('-------------------------');
    disp( seQ0 )
 

    disp('Standard errors (QMLE - p=4)');
    disp('-------------------------');
    disp( seQp )

    
    
    
end


% ======================================================================
% 
% This function calls the lnlt function and returns the average log-likelihood.
%
% ======================================================================

function logl = lnl( b,yt,xt )

    logl = -mean( lnlt( b,yt,xt ) ); 
%    logl = -mean( lnltp( b,y ) ); 
    
end





% ======================================================================
%
%   Define the log of the likelihood at each observation            
% ======================================================================
function loglt = lnlt( b,yt,xt )

    ut    = yt - b(1) - b(2)*xt(:,2) - b(3)*xt(:,3) - b(4)*xt(:,4);
    sig2  = abs( b(5) );
    zt    = ut./sqrt( sig2 );
    loglt = - 0.5*log( 2*pi ) - 0.5*log( sig2 ) - 0.5*zt.^2;


end


% ======================================================================
%
%   Define the log of the likelihood at each observation   
%   Extended for AR(1) errors
%
% ======================================================================
function loglt = lnltp( b,yt,xt )

     ut    = yt - b(1) - b(2)*xt(:,2) - b(3)*xt(:,3) - b(4)*xt(:,4);
     vt    = trimr(ut,1,0) - b(5)*trimr(ut,0,1);
     sig2  = abs( b(6) );
     zt    = vt./sqrt( sig2 );
     loglt = - 0.5*log( 2*pi ) - 0.5*log( sig2 ) - 0.5*zt.^2;

end


