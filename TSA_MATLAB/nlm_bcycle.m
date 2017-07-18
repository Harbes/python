%=========================================================================
%
%   Estimate Markov Switching model (Hamilton, Econometrica, 1989, 357-384)
%
%=========================================================================
function nlm_bcycle( )

    clear all;
    clc;

    % Load data 
    load('nlm_GNP.mat');

    % Compute percentage growth rate (June 1951 to December 1984)	 
    y = 100*( trimr(log(gnp),1,0) - trimr(log(gnp),0,1) );  
    t = length(y);
  
    %**********************************************************************
    %***
    %***     Generate graph
    %***
    %**********************************************************************

    xData = seqa(1951+2/4,1/4,t);
    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;

    plot(xData,y,'-k',xData,1.176*ones(t,1),'--k',xData,-0.224*ones(t,1),'--k');
    ylabel('$y_t$');
    xlabel('t')
    axis([1951 1985 -2 3]) 
    set(gca,'YTick',[-2.0 -1.0 0.0 1.0 2.0 3.0]);
    set(gca,'YTickLabel', {'-2.0','-1.0','0.0','1.0','2.0','3.0'});
    box off

    % Print the tex file to the relevant directory
    %laprint(1,'usgnp','options','factory');

    % Estimate the model without time-varying parameters
    alpha  = 0.0;
    beta   = 0.0;
    gam    = 0.0;
    delta  = 0.0;
    p      = 0.0;
    q      = 0.0;
    % Parameters to control time-variation in the conditional probabilities    
    kappa1 = 0.0;       
    lam1   = 0.0;

    start = [ -0.5, 1, 1.5, 1, 2, 2];
    ops   = optimset('LargeScale','off','Display','off');
    [thetac,lf0] = fminsearch(@(b) neglogr(b,y,1),start,ops);

    clear start
    lf0=-lf0;
    
    % Estimate model again to get standard errros: these are based on the HESSIAN
    start = thetac;   
    start(5) = 1/(1 + exp(-thetac(5)));          
    start(6) = 1/(1 + exp(-thetac(6)));

    [theta,aa,aa,aa,aa,hess] = fminunc(@(b) neglogr(b,y,0),start,ops);
    
    vc = (1/t)*inv(hess);
    se = sqrt(diag(vc));
    
    disp('Markov Switching Model Parameters')
    list  = ['alpha  ';'beta   ';'gamma  ';'delta  ';'p      ';'q      '];
    disp([list, num2str([theta' se])]);

    clear start

    disp(' ')
    disp(['Duration estimate (state=1)     = ',num2str(1/(1 - theta(5)))]);
    disp(['Standard error                  = ',num2str(sqrt(vc(5,5))/(1 - p)^2)]);

    disp(' ')
    disp(['Duration estimate (state=0)     = ',num2str(1/(1 - theta(6)))]);
    disp(['Standard error                  = ',num2str(sqrt(vc(6,6))/(1 - p)^2)]);

    % Wald test (delta = 0)
    w = theta(4)^2/vc(4,4);
    disp(' ')
    disp(['Wald statistic (delta=0)        = ',num2str(w)]);
    disp(['p-value                         = ',num2str(1-chi2cdf(w,1))]);

    % Estimate the model with time-varying probabilities
    start = [ thetac'; 0; 0];
    [theta1,lf1] = fminsearch(@(b) neglog(b,y,urate),start,ops);

    lf1 = -lf1;
    theta1(5) = 1/(1 + exp(-theta1(5)));          
    theta1(6) = 1/(1 + exp(-theta1(6)));

    disp('Time-varying Markov Switching Model Parameters')
    list  = ['alpha  ';'beta   ';'gamma  ';'delta  ';'p      ';'q      ';'kappa  ';'lambda '];
    disp([list, num2str(theta1)]);

    % LR test  
    lr = -2*t*(lf0 - lf1);
    disp(' ')
    disp(['LR statistic (kappa/lambda=0)   = ',num2str(lr)]);
    disp(['p-value                         = ',num2str(1-chi2cdf(lr,2))]);

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Restricted log-likelihood function
%-------------------------------------------------------------------------
function lf = neglogr(b,y,flag)

    % Parameters   
    alpha  = b(1);
    beta   = b(2);
    gam    = b(3);
    delta  = b(4);
    if flag;
     p     = 1/(1 + exp(-b(5)));  % Ensure p,q [0,1]         
     q     = 1/(1 + exp(-b(6)));
    else
        p = b(5);
        q = b(6);
    end

    % Mean and variance vectors for states 1 and 0    
    m  = [alpha + beta  ; alpha ];     
    s2 = [gam   + delta ; gam   ];      

    % Standardized variables for states 1 and 0
    z1 =  (y - m(1))/sqrt(s2(1));      
    z0 =  (y - m(2))/sqrt(s2(2));      
    
    % Distributions for states 1 and 0
    f1 = normpdf(z1)/sqrt(s2(1));         
    f0 = normpdf(z0)/sqrt(s2(2));     

    % Define transitional and stationary probabilities    
     p_stat  = [ (1 - q)/(2 - p - q) ;    
                 (1 - p)/(2 - p - q) ];
     w = p_stat;                        

     % Construct likelihood  
     t = length(y);
     f = zeros(t,1);

     for i = 1:t
        
        % Conditional distribution of y given lagged y 
        f(i) = w(1)*f1(i) + w(2)*f0(i);	
        
        if flag
            p = 1/(1 + exp(-b(5)));           
            q = 1/(1 + exp(-b(6)));
        else
            p = b(5);
            q = b(6);
        end

        p_trans = [   p     1 - q ;
                    1 - p     q   ] ;

        % Update weights of being in each state	
        tmp = [ w(1)*f1(i)/f(i); w(2)*f0(i)/f(i) ];
        w   = p_trans*tmp;             

     end
     
     lf = -mean(log(f) );
end

%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Unrestricted log-likelihood function
%-------------------------------------------------------------------------
function lf = neglog(b,y,urate)

    % Parameters   
    alpha  = b(1);
    beta   = b(2);
    gam    = b(3);
    delta  = b(4);
    p      = 1/(1 + exp(-b(5)));         
    q      = 1/(1 + exp(-b(6)));    
    kappa1 = b(7);   %  cond. prob of being in state 1  
    lam1   = b(8);   %  cond. prob of being in state 0  

    % Mean and variance vectors for states 1 and 0    
    m  = [alpha + beta  ; alpha ];     
    s2 = [gam   + delta ; gam   ];      

    % Standardized variables for states 1 and 0
    z1 =  (y - m(1))/sqrt(s2(1));      
    z0 =  (y - m(2))/sqrt(s2(2));      
    
    % Distributions for states 1 and 0
    f1 = normpdf(z1)/sqrt(s2(1));         
    f0 = normpdf(z0)/sqrt(s2(2));     

    % Define transitional and stationary probabilities    
    p_stat  = [ (1 - q)/(2 - p - q) ;    
                 (1 - p)/(2 - p - q) ];
    w       = p_stat;                        
    
    % Construct likelihood  
    t = length(y);
    f = zeros(t,1);
        
    for i = 1:t
        
        % Conditional distribution of y given lagged y 
        f(i) = w(1)*f1(i) + w(2)*f0(i);		                                           
        p    = 1/(1 + exp(-b(5)-kappa1*urate(i)));           
        q    = 1/(1 + exp(-b(6)-lam1*urate(i)));

        p_trans = [   p     1 - q ;
                    1 - p     q   ] ;

        % Update weights of being in each state	
        tmp = [ w(1)*f1(i)/f(i); w(2)*f0(i)/f(i) ];
        w   = p_trans*tmp;             

    end

    lf = -mean(log(f) ) ;

end




