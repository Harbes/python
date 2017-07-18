%========================================================================
%
%      Estimate investment model and compute QMLE standard errors
%
%========================================================================

function qmle_invest_auto( )

    clear all
    clc

    % Load quarterly data for the US for the period 1957 to 2007
    load usinvest.mat

    % Generate variables: data start in 1958Q1 
    ri = log(trimr(invest./cpi,1,0)); 
    ry = log(trimr(gdp./cpi,1,0));                    
    inf = log(trimr(cpi,1,0)./trimr(cpi,0,1));            
    rint = trimr(r10yr/100,1,0) - 4*inf;                            

    t  = length(ri);
    yt = ri;
    xt = [ ones(t,1) ry rint ];

    
    % Estimate model by OLS
    bols = xt\yt;
    et   = yt - xt*bols;
    s2   = mean(et.^2);

    
    % Find parameters using the function with rho1 constrained
    bstart = [ bols; tanh(0.998); s2 ];
    options = optimset( 'LargeScale','off',  ...
                        'Display','iter',...
                        'MaxFunEvals',5000, ...
                        'MaxIter',2000);   

    [ bhat ] = fminunc( @(b) lnlc( b,yt,xt ),bstart,options);
    
    bhat(4) = atanh(bhat(4));
    
    % Estimate the gradient and Hessian using the unconstrained parameters    
    G = numgrad( @lnlt,bhat,yt,xt );
    H = numhess( @lnl,bhat,yt,xt );
    
    % Standard Errors based on Hessian
    iH  = inv(H);
    seH = sqrt( diag( iH ) );
    
   
    % Compute standard errors based on the OPG matrix
    J   = G'*G;
    seG = sqrt( diag( inv( J ) ) );
 
    % Compute QMLE standard errors  
    seQ0 = sqrt( diag( iH*J*iH) );

 
    %  Compute qmle standard errors with lag length based on p 
    p = floor( 4*(t/100)^(2/9) );
 
    t = size(G,1);
    for i = 1:p;
 
          gmat = G((i+1):t,:)'*G(1:(t-i),:);
          J    = J + (1.0 - i/(p+1))*(gmat + gmat');
    
    end

    
    seQp = sqrt( diag( iH*J*iH) );
     

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

%
%--------------------------- Functions -----------------------------------
% 

%-------------------------------------------------------------------------
% Log-likelihood function at each observation with constraint on rho1
%-------------------------------------------------------------------------
function logl = lnlc( b,yt,xt )

    ut   = yt - b(1) - b(2)*xt(:,2) - b(3)*xt(:,3);
    vt   = trimr(ut,1,0) - atanh(b(4))*trimr(ut,0,1);
    sig2 = abs( b(5) );
    zt   = vt./sqrt( sig2 );
    tmp  = -0.5*log( 2*pi ) - 0.5*log( sig2 ) - 0.5*zt.^2;
    logl = -sum( tmp );

end

%-------------------------------------------------------------------------
% Wrapper function for unconstrained log-likelihood
%-------------------------------------------------------------------------


function logl = lnl( b,yt,xt )

    logl = -sum( lnlt( b,yt,xt ) ); 
    
end


% ======================================================================
%
%   Define the log of the likelihood at each observation   
%   Extended for AR(1) errors
%
% ======================================================================
function logl = lnlt( b,yt,xt )

    ut   = yt - b(1) - b(2)*xt(:,2) - b(3)*xt(:,3);
    vt   = trimr(ut,1,0) - b(4)*trimr(ut,0,1);
    sig2 = b(5);
    zt   = vt./sqrt( sig2 );
    logl = -0.5*log( 2*pi ) - 0.5*log( sig2 ) - 0.5*zt.^2;

end


