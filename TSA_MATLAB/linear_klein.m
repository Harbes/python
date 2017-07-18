% ========================================================================
% 
%   Program to estimate Klein's macroeconomic model by full 
%   information maximum likelihood.
% 
% ========================================================================
function linear_klein(  ) 		 

    clear all;
    clc;
    
    t = 22;

    % Read in the data   
    klein_data = load('klein.dat', '-ASCII' );  
    
    c      = klein_data(:,1);
    p      = klein_data(:,2);
    pw     = klein_data(:,3);
    i      = klein_data(:,4);
    klag   = klein_data(:,5);
    d      = klein_data(:,6);
    income = klein_data(:,7);
    gw     = klein_data(:,8);
    g      = klein_data(:,9);
    tax    = klein_data(:,10);
    trend  = (-11:1:-11+t-1)';

    % Estimate consumption function by OLS        
    y         = c(2:end,:);
    pwgw      = pw + gw;
    x         = [ones(t-1,1) p(2:end,:) p(1:end-1) pwgw(2:end,:)];
    alpha_ols = x\y;

    disp('OLS parameter estimates of the consumption function:');
    disp( alpha_ols );

    % Estimate investment function by OLS        
    y        = i(2:end,:);
    x        = [ones(t-1,1) p(2:end,:) p(1:end-1) klag(2:end,:)];
    beta_ols = x\y;
    
    disp('OLS parameter estimates of the investment function:');
    disp( beta_ols );

    % Estimate wage function by OLS        
    y       = pw(2:end,:);
    x       = [ones(t-1,1) d(2:end,:) d(1:end-1) trend(2:end,:)];
    gam_ols = ((x'*y)'/(x'*x));

    disp('OLS parameter estimates of the wage function:');
    disp( gam_ols );

    % Define the endogenous variables for the system  
    cipw = [c i pw];
    y    = cipw(2:end,:);

    % Define exogenous/predetermined variables (instruments) for the system    
    x = [ones(t-1,1) g(2:end,:) tax(2:end,:) gw(2:end,:) trend(2:end,:) p(1:end-1,:) d(1:end-1,:) klag(2:end,:)];

    % Estimate consumption function by IV        
    pwgw = pw + gw;
    alpha_iv = iv(c(2:end,:), [ones(t-1,1) p(2:end,:) p(1:end-1,:) pwgw(2:end,:)], x);

    disp('IV estimates of the consumption function');
    disp( alpha_iv );

    % Estimate investment function by IV       
    beta_iv = iv(i(2:end,:), [ones(t-1,1) p(2:end,:) p(1:end-1,:) klag(2:end,:)], x);

    disp('IV estimates of the investment function');
    disp( beta_iv );

    % Estimate wage function by IV        
    gam_iv = iv(pw(2:end,:), [ones(t-1,1) d(2:end,:) d(1:end-1) trend(2:end,:)], x);
    
    disp('IV estimates of the wage function');
    disp( gam_iv );

    % Estimate the model by FIML using IV starting values
    theta0 = [alpha_iv; beta_iv; gam_iv;];       
    [theta,a0]  = fminunc(@(theta) neglog(theta,y,x),theta0);

    disp('FIML estimates of the consumption function');
    disp( theta(1:4) );
    disp('FIML estimates of the investment function');
    disp( theta(5:8) );
    disp('FIML estimates of the wage function');
    disp( theta(9:12) );
    disp( ' ' );
    disp( ['Log-likelihood = ', num2str(-a0) ] );

    [t n] = size(y);
    b = [1-theta(2)		       -theta(6)	-theta(10);
          -theta(2) 		  1-theta(6) 	-theta(10);
           theta(2)-theta(4)   -theta(6) 		1];

    a = [-theta(1) 		-theta(5) 		-theta(9);
         -theta(2) 		-theta(6) 		-theta(10);
          theta(2) 		theta(6) 		0;
         -theta(4) 		0 			    0;
         0 				0 			   -theta(12);
        -theta(3) 		-theta(7) 		0;
         0 				0 			   -theta(11);
         0 				-theta(8) 		0];

     % Compute residuals and covariance matrix
    u = zeros(t,n);
    for i=1:t
        u(i,:) = y(i,:)*b + x(i,:)*a;     	                  
    end
    rvc  = u'*u/t;                       	

    disp('Estimated covariance matrix');
    disp( rvc );

end
%
%--------------------------- Functions -----------------------------------
% 

%-------------------------------------------------------------------------
% Log-likelihood function 
%-------------------------------------------------------------------------
function lf = neglog( theta,y,x )

    lf = -mean( lnlt( theta,y,x ) );

end
%-------------------------------------------------------------------------
% Log-likelihood function at each observation     
%-------------------------------------------------------------------------
function lnl = lnlt(theta,y,x);
	
    [t n] = size(y);
    
    b = [1-theta(2)		       -theta(6)	-theta(10);
          -theta(2) 		  1-theta(6) 	-theta(10);
           theta(2)-theta(4)   -theta(6) 		1];

    a = [-theta(1) 		-theta(5) 		-theta(9);
         -theta(2) 		-theta(6) 		-theta(10);
          theta(2) 		theta(6) 		0;
         -theta(4) 		0 			    0;
         0 				0 			   -theta(12);
        -theta(3) 		-theta(7) 		0;
         0 				0 			   -theta(11);
         0 				-theta(8) 		0];
         
    % Construct residuals and covariance matrix
    u = zeros(t,n);
    for i=1:t
        u(i,:) = y(i,:)*b + x(i,:)*a;     	                    
    end
    omega = u'*u/t;                       	  
    
    lnl = zeros(t,1);    
   	for i=1:t
	    lnl(i) = -n*0.5*log(2*pi) + log(abs(det(b))) - 0.5*log(det(omega)) ...
            - 0.5*u(i,:)*inv(omega)*u(i,:)';
    end
    
end	
	
%-------------------------------------------------------------------------
% Instrumental variable estimation
%-------------------------------------------------------------------------
function b = iv(y,w,x)

    % IV estimates
    b      = inv(w'*x*inv(x'*x)*x'*w)*(w'*x*inv(x'*x)*x'*y);    	                            
	
    % % Standard error of regression
    e      = y - w*b;                                           	                          	
    k      = size(w,2);                                           	           
	t      = size(y,1);
	sigma  = sqrt( e'*e/t );                                     	            
	
    % Variance-covariance matrix
    vcov   = sigma^2*inv(w'*x*inv(x'*x)*x'*w);                     	            
	sterr  = sqrt( diag(vcov) );                                	                        
	tstats = b./sterr;                                          	                          
    
end		