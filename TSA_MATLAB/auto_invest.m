%=========================================================================
%
%     Program to estimate a dynamic investment model
%     U.S. quarterly data for period 1957 to 2007
%
%=========================================================================

function auto_invest( )

    clear all;
    clc;

    % Load data
    load usinvest.txt
    
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
     
    % Exact MLE
    theta0  = [b ; 0.02  ; sqrt(s2)] ;         
    [theta,fe,~,~,~,H]  = fminunc(@(p) neglog(p,dri,dry,rint),theta0);
    
    invH = inv(H);
     
    disp('Results for exact MLE');
    disp(['Log-likelihood function = ', num2str(-fe)]);
    disp('Parameter estimates and std. errors');
    disp( [theta diag(invH)/t] );
  
    % Conditional MLE 
    [theta,fc,~,~,~,H]  = fminunc(@(p) neglogc(p,dri,dry,rint),theta);

    invH = inv(H);
    
    disp('Results for conditional MLE');
    disp(['Log-likelihood function = ', num2str(-fc)]);
    disp('Parameter estimates and std. errors');
    disp( [theta diag(invH)/t] );
    
end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%   Exact log-likelihood function 
%-------------------------------------------------------------------------
function lf = neglog(b,ri,ry,rint)

     beta0 = b(1);                                          
     beta1 = b(2);
     beta2 = b(3);
     rho1  = b(4);
     sig2  = b(5);

     u     = ri - beta0 - beta1*ry - beta2*rint;
     v     = trimr(u,1,0) - rho1*trimr(u,0,1);
     
     %  Log-likelihood for t=1              
     lnl_0 = - 0.5*log(2*pi) - 0.5*log(sig2) + 0.5*log(1 - rho1^2) ...
         - 0.5*(u(1) - 0).^2/(sig2/(1 - rho1^2));   
     %  Log-likelihood for t=2,3,...   
     lnl_1 = - 0.5*log(2*pi) - 0.5*log(sig2) - 0.5*v.^2/sig2;                                                  

     lf = -mean( [lnl_0 ; lnl_1] );
     
end
%-------------------------------------------------------------------------
%  Conditional log-likelihood function 
%-------------------------------------------------------------------------
function lf = neglogc(b,ri,ry,rint)

     beta0 = b(1);                                          
     beta1 = b(2);
     beta2 = b(3);    
     rho1  = b(4);
     sig2  = b(5);

     u     = ri - beta0 - beta1*ry - beta2*rint;
     v     = trimr(u,1,0) - rho1*trimr(u,0,1);
     
     %  Log-likelihood for t=2,3,...   
     lnl_1 = - 0.5*log(2*pi) - 0.5*log(sig2) - 0.5*v.^2/sig2;                                                  

     lf = -mean( lnl_1 );
     
end

    