%=========================================================================
%
%     Program to estimate and test a dynamic investment model
%     U.S. quarterly data for period 1957 to 2007
%     
%=========================================================================

function auto_test( )

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
  
    % Unconstrained model
    theta  = [b ; 0.02  ; sqrt(s2)] ;         
    [theta1,f1,~,~,~,H1]  = fminunc(@(p) neglog1(p,dri,dry,rint),theta);

    lnl1  = -f1;
    invH1 = inv(H1);
    
    disp('Results for conditional MLE');
    disp(['Log-likelihood function = ', num2str(lnl1)]);
    disp('Parameter estimates and std. errors');
    disp( [theta diag(invH1)/t] );
    
    % Constrained model 
    theta = [b ; sqrt(s2)] ;
    [theta0,f0,~,~,~,H0]  = fminunc(@(p) neglog0(p,dri,dry,rint),theta);
    
    invH0 = inv(H0);
    lnl0  = -f0;
    
    disp('Results for conditional MLE');
    disp(['Log-likelihood function = ', num2str(lnl0)]);
    disp('Parameter estimates and std. errors');
    disp( [theta diag(invH0)/t] );

    % LR test      
    lr = -2*t*(lnl0 - lnl1);
    disp(['LR statistic            = ',num2str(lr) ]); 
    disp(['p-value                 = ',num2str(1-cdf('chi2',lr,1)) ]);

    % Wald test   
    r = [ 0 , 0 , 0 , 1 , 0 ];
    q = 0;
    wd = t*(r*theta1 - q)'*inv(r*inv(H1)*r')*(r*theta1 - q);
    disp(['Wald statistic          = ',num2str(wd) ]); 
    disp(['p-value                 = ',num2str(1-cdf('chi2',wd,1)) ]);

    % LM test      
    theta = [ theta0(1:3) ; 0.0 ; theta0(4) ];
    gmat  = numgrad(@lnlt1,theta,dri,dry,rint);  
    g     = mean(gmat)';
    j     = gmat'*gmat/t;
    lm    = t*g'*inv(j)*g; 
    disp('Gradient evaluated at contrained estimates');
    disp( g );
	disp('Outer product of gradients matrix');
    disp( j );     
    disp(['LM statistic            = ',num2str(lm) ]); 
    disp(['p-value                 = ',num2str(1-cdf('chi2',lm,1)) ]);                

    % LM test (regression form)                                
    % Stage 1 regression
    b = x\y; 
    u = y - x*b;                                  
    
    % Stage 2 regression
    y = trimr(u,1,0);
    z = [trimr(x,1,0) , trimr(e,0,1)];
    
    b = z\y; 
    e  = y - z*b;
    r2= 1 - ((y-mean(y))'*(y-mean(y)))\e'*e;
    lm = (t-1)*r2;
    disp(['LM statistic (regression) = ',num2str(lm) ]); 
    disp(['p-value                   = ',num2str(1-cdf('chi2',lm,1)) ]);

    %  LM test (first-order autocorrelation of residuals)           
    y = dri;
    b = x\y;                    
    u = y - x*b;                  
    r1 = autocorr(u,1,0);     
    lm = (t-1)*r1(2)^2;

    disp(['LM statistic (alternative) = ',num2str(lm) ]); 
    disp(['p-value                    = ',num2str(1-cdf('chi2',lm,1)) ]);

end
%
%--------------------------- Functions -----------------------------------
% 
%-----------------------------------------------------------------------
% Negative unconstrained log-likelihood  
%-----------------------------------------------------------------------
function lf = neglog1(b,ri,ry,rint)

    lf = -mean( lnlt1(b,ri,ry,rint) );

end
%-----------------------------------------------------------------------
% Unconstrained log-likelihood function at each observation
%-----------------------------------------------------------------------
function lnl = lnlt1(b,ri,ry,rint)

     beta0 = b(1);                                          
     beta1 = b(2);
     beta2 = b(3);    
     rho1  = b(4);
     sig2  = b(5);

     u     = ri - beta0 - beta1*ry - beta2*rint;
     v     = trimr(u,1,0) - rho1*trimr(u,0,1);
     
     %  Log-likelihood for t=2,3,...   
     lnl = - 0.5*log(2*pi) - 0.5*log(sig2) - 0.5*v.^2/sig2;                                                  
     
end
%-----------------------------------------------------------------------
% Negative constrained log-likelihood function
%-----------------------------------------------------------------------
function lf = neglog0(b,ri,ry,rint)

     beta0 = b(1);                                                     
     beta1 = b(2);
     beta2 = b(3);
     rho1  = 0.0;
     sig2  = b(4);
     u     = ri - beta0 - beta1*ry - beta2*rint;
     v     = trimr(u,1,0) - rho1*trimr(u,0,1);
     
     % Log-likelihood for t=2,3,...
     lnl = - 0.5*log(2*pi) - 0.5*log(sig2) - 0.5*v.^2/sig2;     
     lf  = - mean( lnl );
end
     