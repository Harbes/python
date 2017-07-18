% =======================================================================
%
%   Program to estimate a model by full information maximum likelihood 
% 	and instrumental variables. The set of equations is given by the 
%   bivariate system
%
%         y1t = beta*y2t + u1t
%         y2t = gam*y1t + alpha*xt + u2t
%
% =======================================================================

function linear_iv(  ) 		

    clear all;
    clc;

    % Get data
    t       = 500;   
    [ y,x ] = simulatedata( t );

    % Estimate the model using random starting values
    theta  = fminunc( @(theta) neglog(theta,y,x),rand(3,1) );
    
    disp(' MLE Estimates ' );
    disp( theta );


    % Construct residuals and estimate covariance matrix
    [t n] = size(y);
    b = [1 -theta(2); -theta(1) 1]; 
    a = [0 -theta(3)];
    e = zeros(t,n);

    for i=1:t
        
        e(i,:) = y(i,:)*b + x(i,:)*a;     
        
    end
    me2   = mean(e.^2);
    drv   = eye(n);
    omega = diag( me2 );
    
    disp('Covariance matrix');
    disp( omega );
      
    % FIML analytical solution    

    beta_fiml  = sum(y(:,1).*x) / sum(y(:,2).*x);
    e1         = y(:,1) - beta_fiml*y(:,2);
    num        = sum(y(:,2).*e1)*sum(x.^2) - sum(x.*e1)*sum(y(:,2).*x);
    den        = sum(y(:,1).*e1)*sum(x.^2) - sum(x.*e1)*sum(y(:,1).*x);
    gam_fiml   = num/den;
    num        = sum(y(:,1).*e1)*sum(y(:,2).*x) - sum(y(:,1).*x)*sum(y(:,2).*e1);
    alpha_fiml = num/den;

    disp('Closed form FIML estimates') 
    disp( [beta_fiml; gam_fiml; alpha_fiml ] );



    % Instrumental variable estimation    
    beta1_iv = iv(y(:,1), y(:,2), x);       
    beta2_iv = iv(y(:,2), [y(:,1) x], [e1 x]);

    disp('Instrumental variables estimates')
    disp( [beta1_iv; beta2_iv(1); beta2_iv(2) ] );

end

%
%--------------------- Functions ----------------------%
%
% Simulate data
function [ y,x ] = simulatedata( t )

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123457) );

    beta  = 0.6; 
    gam   = 0.4;
    alpha = -0.5;  
    omega = [2.0 0.0; 0.0 1.0];

    % Construct population parameter matrices     
    b = [1 -gam; -beta 1]; 
    a = [0 -alpha];

    % Construct exogenous variables                       
    x = 10*randn(t,1);

    % Construct disturbances                              
    u = randn(t,2)*chol(omega);

    % Simulate the model by simulating the reduced form   
    y = zeros(t,2);

    for i=1:t
        y(i,:) = -x(i,:)*a*inv(b) + u(i,:)*inv(b);
    end
end

% Negative log-likelihood function
function lf = neglog( theta,y,x )

    lf = -mean( lnlt(theta,y,x) );

end

% Log-likelihood at each observation         
function lnl = lnlt(theta,y,x)

    [t n] = size(y);
    
    b = [1 -theta(2); -theta(1) 1]; 
    a = [0 -theta(3)];
    e = zeros(t,n);

    % Construct residuals
    for i=1:t
        
        e(i,:) = y(i,:)*b + x(i,:)*a;     	                    
    
    end

	% Concentrate out resid var-covar matrix and restrict it to be diagonal	
	me2 = mean(e.^2);
	omega = eye(n);
    omega = diag( me2 );
    
    lnl = zeros(t,1);
    
   	for i=1:t
	    lnl(i) = -n*0.5*log(2*pi) + log(abs(det(b))) - 0.5*log(det(omega)) ...
            - 0.5*e(i,:)*inv(omega)*e(i,:)';
    end
    
end	

% Instrumental variable estimation
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