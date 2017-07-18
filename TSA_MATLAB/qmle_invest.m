%========================================================================
%
%      Estimate investment model and compute QMLE standard errors
%
%========================================================================

function qmle_invest( )

    clear all;
    clc;

    % Load data
    load investment
    
    cpi    = usinvest(:,1); 
    gdp    = usinvest(:,2);
    invest = usinvest(:,3);
    r10yr  = usinvest(:,4);
    r3yr   = usinvest(:,5);
    tbill  = usinvest(:,6);
     
    % Generate variables: data start in 1958Q1 
    
    gfc  = [zeros(length(usinvest)-13,1) ; ones(13,1) ] ; 
    dri  = 100*(log(trimr(invest./cpi,1,0)) - log(trimr(invest./cpi,0,1)));   
    inf  = 100*log(trimr(cpi,1,0)./trimr(cpi,0,1));                           
    rint = trimr(r10yr/4,1,0) - inf;                                          
    dry  = 100*( log(trimr(gdp./cpi,1,0)) - log(trimr(gdp./cpi,0,1)) );       
    gfc  = trimr(gfc,1,0);
    
    % OLS regression
    t = length(dry);
    y = dri;
    x = [ones(t,1),dry,rint];

    b  = x\y;
    e  = y - x*b;
    s2 = mean(e.^2);
      
    % Conditional MLE (same as the OLS estimates)
    theta0 = [ b  ; s2 ];
    [theta,fc,~,~,~,H]  = fminunc(@(p) neglog(p,dri,dry,rint),theta0);

    iH  = inv(H);
    seH = (1/t)*diag(iH);
    
    disp('Results for conditional MLE');
    disp(['Log-likelihood function = ', num2str(-fc)]);
    disp('Parameter estimates and std. errors');
    disp( [theta seH] );
           
    % Estimate the gradient and Hessian using the unconstrained parameters 
    g    = numgrad( @lnlt,theta,dri,dry,rint );
    
   
    % Compute standard errors based on the OPG matrix
    j   = g'*g/t;
    seG = sqrt( (1/t)*diag( inv( j ) ) );
 
    % Compute QMLE standard errors  
    seQ0 = sqrt( (1/t)*diag( iH*j*iH) );

 
    %  Compute qmle standard errors with lag length based on p 
    p = floor( 4*(t/100)^(2/9) );
 
    tmax = size(g,1); 
    for i = 1:p;
 
          gmat = g((i+1):tmax,:)'*g(1:(tmax-i),:)./t;
          j    = j + (1.0 - i/(p+1))*(gmat + gmat');
    
    end

    
    seQp = sqrt( (1/t)*diag( iH*j*iH) );
     

    disp('Parameter Estimates');
    disp('-------------------');
    disp( theta );
 
    
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
% Negative unconstrained log-likelihood function
%-------------------------------------------------------------------------


function logl = neglog( b,ri,ry,rint )

    logl = -mean( lnlt( b,ri,ry,rint ) ); 
    
end


%-------------------------------------------------------------------------
%  Conditional log-likelihood function 
%-------------------------------------------------------------------------
function lf = lnlt(b,ri,ry,rint)

     beta0 = b(1);                                          
     beta1 = b(2);
     beta2 = b(3);    
     sig2  = b(4);

     u     = ri - beta0 - beta1*ry - beta2*rint;
     
     %  Log-likelihood for t=2,3,...   
     lf = - 0.5*log(2*pi) - 0.5*log(sig2) - 0.5*u.^2/sig2;                                                  
    
end



