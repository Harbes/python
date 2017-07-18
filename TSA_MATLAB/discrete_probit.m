%=========================================================================
%
%   Program to estimate a probit model of US monetary policy
%   The data file is based on the Hamilton and Jorda (2002) data set.
%
%=========================================================================
function discrete_probit( )

    clear all
    clc

    % Read the data: US weekly yields   
    % from 1st week of February 1984 to first week of June 1997
    load usmoney.mat
    event   = usmoney(:,1);
    target  = usmoney(:,2);
    change  = usmoney(:,3);
    bin     = usmoney(:,4);
    fomc    = usmoney(:,8);
    spread6 = usmoney(:,20);
    inf     = usmoney(:,23);
    unemp   = usmoney(:,25);
    gdp     = usmoney(:,28);
    
    % Reverse the spread so it is the Federal funds rate 
    % less 6-month Treasury bill rate       
    spread = -spread6;                                                  

    %  Redefine the target rate based on the consolidated series 
    % constructed in Hamilton and Jorda (2002)       
    target_adj = cumsum( [ target(1) ; bin(2:end) ] );

    plot(seqa(1984+5/52,1/52,length(target_adj)),target_adj);
    
    % Choose data based on fomc days
    ind       = fomc == 1;
    data      = [ bin spread inf gdp ];
    data_fomc = data(ind,:);
    
    % Dependent and independent variables 
    y = double( data_fomc(:,1) > 0.0 );
    t = length(y);
    x = [ ones(t,1) data_fomc(:,2:4) ];   
    
    % Estimate model by OLS (ie ignore that the data are binary
    b = x\y;
    u = y - x*b;
    s = sqrt( mean(u.^2) );

    % Estimate the unrestricted probit regression model by MLE   
    ops    = optimset('LargeScale','off','Display','iter');
    theta0 = b/s;

    [theta1,l1,~,~,~,h] = fminunc(@(b) lprobit(b,y,x),theta0,ops);
 
    cov = (1/t)*inv(h);
    l1 = -l1;
    
    disp(['Unrestricted log-likelihood function =     ',num2str(l1) ]);
    disp(['T x unrestricted log-likelihood function = ',num2str(t*l1)]);
    
    disp('Unrestricted parameter estimates')
    disp(theta1)
 
    % Estimate the restricted probit regression model by MLE     
    theta0 = norminv(mean(y),0,1);
    [theta,l0,~,~,~,h] = fminunc(@(b) l0probit(b,y),theta0,ops);
    
    l0 = -l0;
    disp(['Unrestricted log-likelihood function =     ',num2str(-l0) ]);
    disp(['T x unrestricted log-likelihood function = ',num2str(-t*l0)]);

    disp('Restricted parameter estimates')
    disp(theta)
 

    %  Likelihood ratio test  
	lr = -2*t*(l0 - l1);
 
	disp(['LR Statistic         = ',num2str(lr)]);
    disp(['p-value              = ',num2str(1-chi2cdf(lr,size(x,2)-1))]);
    disp(' ');
    
    % Wald test  
    r = [ 0   1   0   0 ;
          0   0   1   0 ;
          0   0   0   1 ];

    q = [ 0 ; 0 ; 0 ];

    wd = (r*theta1 - q)'*inv(r*cov*r')*(r*theta1 - q);

	disp(['Wald Statistic       = ',num2str(wd)]);
    disp(['p-value              = ',num2str(1-chi2cdf(wd,size(x,2)-1))]);

    disp(' ');


    % LM test of the joint restrictions  
    u  = y - mean(y);
    b  = x\u;
    e  = u - x*b;
    r2 = 1 - (e'*e)/(u'*u);
    lm = t*r2;

    disp(['Regression estimates = ',num2str( b') ]);
    disp(['Sample size (t)      = ',num2str( t ) ]);
    disp(['R2                   = ',num2str( r2 ) ]);
	disp(['LM Statistic         = ',num2str(lm)]);
    disp(['p-value              = ',num2str(1-chi2cdf(lm,size(x,2)-1))]);
    
end

%
%--------------------------- Functions -----------------------------------
% 

%-------------------------------------------------------------------------
%  Unrestricted Probit negative log-likelihood function 
%-------------------------------------------------------------------------
function lf = lprobit(b,y,x)


    f  = normcdf( x*b ); 
    lf = -mean( y.*log(f) + (1 - y).*log(1 - f) );

end 

%-------------------------------------------------------------------------
%  Restricted Probit negative log-likelihood function 
%-------------------------------------------------------------------------
function lf = l0probit(b,y)


    f  = normcdf( b ); 
    lf = -mean( y.*log(f) + (1 - y).*log(1 - f) );

end 

