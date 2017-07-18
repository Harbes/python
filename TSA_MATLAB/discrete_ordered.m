%=========================================================================
%
%   Program to estimate an ordered probit model of US monetary policy
%   The data file is based on the Hamilton and Jorda (2002) data set.
%
%=========================================================================
function discrete_ordered( )

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
    
    % Choose data based on fomc days
    ind       = fomc == 1;
    data      = [ bin spread inf gdp ];
    data_fomc = data(ind,:);
    
    % Dependent and independent variables 
    y = data_fomc(:,1);
    t = length(y);
    x = [ ones(t,1) data_fomc(:,2:4) ];   
    
    % Create dummy variables for each interest rate change
    d1 = double(y == -0.50);
    d2 = double(y == -0.25);
    d3 = double(y ==  0.00);
    d4 = double(y ==  0.25);
    d5 = double(y ==  0.50);
    
    d  = [ d1 d2 d3 d4 d5];
    
    % Estimate model by OLS (ie ignore that the data are ordered
    b = x\y;
    u = y - x*b;
    s = sqrt( mean(u.^2) );
    
    % Compute the unconditional probabilities of each regime    
    p =  cumsum( mean(d) ); 
    
    % Estimate the unrestricted ordered probit regression model by MLE   
    ops    = optimset('LargeScale','off','Display','iter');
    theta0 = [norminv(p(1:4),0,1)' ; b(2:4)/s ];

    [theta1,l1,~,~,~,h] = fminunc(@(b) lprobit(b,x,d),theta0,ops);
    
    disp('Unrestricted parameter estimates')
    disp( theta1 )

    cov = (1/t)*inv(h);
    l1 = -l1;
    
    disp(['Unrestricted log-likelihood function =     ',num2str(l1) ]);
    disp(['T x unrestricted log-likelihood function = ',num2str(t*l1)]);
 
    % Estimate the restricted probit regression model by MLE     
    theta0 = norminv(p(1:4),0,1)';
    [theta,l0,~,~,~,h] = fminunc(@(b) l0probit(b,d),theta0,ops);
  
    disp('Restricted parameter estimates')
    disp( theta )

    
    l0 = -l0;
    disp(['Unrestricted log-likelihood function =     ',num2str(-l0) ]);
    disp(['T x unrestricted log-likelihood function = ',num2str(-t*l0)]);

    %  Likelihood ratio test  
	lr = -2*t*(l0 - l1);
 
	disp(['LR Statistic         = ',num2str(lr)]);
    disp(['p-value              = ',num2str(1-chi2cdf(lr,size(x,2)-1))]);

    disp(' ');
    
    % Wald test  
    r = [ 0   0   0   0   1   0   0 ;
          0   0   0   0   0   1   0 ;
          0   0   0   0   0   0   1 ];


    q = [ 0 ; 0 ; 0 ];

    wd = (r*theta1 - q)'*inv(r*cov*r')*(r*theta1 - q);

	disp(['Wald Statistic       = ',num2str(wd)]);
    disp(['p-value              = ',num2str(1-chi2cdf(wd,size(x,2)-1))]);

    disp(' ');

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Unrestricted Probit negative log-likelihood function 
%-------------------------------------------------------------------------
function lf = lprobit(b,x,d)

    % Cut off points
    c = b(1:4);
    
    % Regression part excluding the intercepts
    xb = x(:,2:4)*b(5:7);
    
    % Cut off points
    f1 = normcdf( c(1) - xb,0,1);
	f2 = normcdf( c(2) - xb,0,1) - f1;
	f3 = normcdf( c(3) - xb,0,1) - f1 - f2;
	f4 = normcdf( c(4) - xb,0,1) - f1 - f2 - f3;
	f5 = 1                       - f1 - f2 - f3 - f4;
    f  = [ f1  f2  f3  f4  f5 ];   
    
    % Log-likelihood function
    tp = bsxfun(@times,d,log(f));
    lf = -mean( sum(tp,2) );     

end 
%-------------------------------------------------------------------------
%  Restricted Probit negative log-likelihood function 
%-------------------------------------------------------------------------
function lf = l0probit(b,d)
       
    % Cut off points
    f1 = normcdf( b(1),0,1);
	f2 = normcdf( b(2),0,1) - f1;
	f3 = normcdf( b(3),0,1) - f1 - f2;
	f4 = normcdf( b(4),0,1) - f1 - f2 - f3;
	f5 = 1                       - f1 - f2 - f3 - f4;
    f  = [ f1  f2  f3  f4  f5 ];   
    
    % Log-likelihood function
    tp = bsxfun(@times,d,log(f));
    lf = -mean( sum(tp,2) );     

end 


