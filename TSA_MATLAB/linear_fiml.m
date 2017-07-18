%=========================================================================
%
%    Program to estimate model by full information maximum likelihood.
%
%=========================================================================
function linear_fiml( ) 		

    clear all;
    clc;
    t = 500;
    %flag = 1;       % 1 = simulate data
    flag = 0;        % 0 = use GAUSS data to reproduce text
   
    if flag 
        [ y,x ] = simulatedata( t );
    else      
        % Aternatively load the GAUSS data
        load linear_fimldata.dat
        y = linear_fimldata(:,[1 2]);
        x = linear_fimldata(:,[3 4]);     
    end
    
    % Estimate the model with random starting values
    theta   = fminunc(@(theta) neglog(theta,y,x),rand(4,1)); 
    ht      = numhess(@neglog,theta,y,x);
    cov     = (1/t)*inv(ht);
    u       = zeros(t,2);
    for i=1:t
        
        u(i,:) = y(i,:)*[1 -theta(3); -theta(1) 1] + x(i,:)*[-theta(2) 0; 0 -theta(4)];     	                    
    
    end

    disp(['Beta 1 and se = ', num2str(theta(1)),'   ', num2str(sqrt(cov(1,1))) ]);
    disp(['Beta 2 and se = ', num2str(theta(2)),'   ', num2str(sqrt(cov(2,2)))]);
    disp(['Beta 3 and se = ', num2str(theta(3)),'   ', num2str(sqrt(cov(3,3)))]);
    disp(['Beta 4 and se = ', num2str(theta(4)),'   ', num2str(sqrt(cov(4,4)))]);

    disp('Residual variance-covariance matrix')
    disp(u'*u/t);

    % Instrumental variable estimation    
    beta1_iv = iv( y(:,1), [y(:,2) x(:,1)], x(:,[1 2]) );

    disp(['IV estimate of beta1   = ', num2str(beta1_iv(1))]);
    disp(['IV estimate of alpha1  = ', num2str(beta1_iv(2))]);

    iv( y(:,2), [y(:,1) x(:,2)], x(:,[1 2]) );

    beta2_iv = iv( y(:,2), [y(:,1) x(:,2)], x(:,[1 2]) );

    disp(['IV estimate of beta2   = ', num2str(beta2_iv(1))]);
    disp(['IV estimate of alpha2  = ', num2str(beta2_iv(2))]);

end
%
%--------------------- Functions ----------------------%
%
% Simulate data
function [ y,x ] = simulatedata( t )

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123457) );

    beta1  = 0.6; 
    alpha1 = 0.4;
    beta2  = 0.2; 
    alpha2 = -0.5;  
    omega  =  [1 0.5; 0.5 1.0];

    % Construct population parameter matrices     
    b = [1 -beta2; -beta1 1]; 
    a = [-alpha1 0; 0 -alpha2];

    % Construct exogenous variables                       
    x = [10*randn(t,1) 3*randn(t,1)];

    % Construct disturbances                              
    u = randn(t,2)*chol(omega);

    % Simulate the model by simulating the reduced form   
    y = zeros(t,2);

    for i=1:t
        y(i,:) = -x(i,:)*a*inv(b) + u(i,:)*inv(b);
    end
end

% Log-likelihood function            
function lf = neglog(theta,y,x)
    
    lf = -mean(lnlt(theta,y,x));
    
end

% Log-likelihood function at each observation            
function rval = lnlt(theta,y,x)
	
    
    [t n] = size(y);
    
    b = [1 -theta(3); -theta(1) 1]; 
    a = [-theta(2) 0; 0 -theta(4)];
    u = zeros(t,n);

    for i=1:t
        u(i,:) = y(i,:)*b + x(i,:)*a;     	                    
    end

    omega = u'*u/t;                         % Concentrate out resid var-covar matrix  
    lnl = zeros(t,1);

    for i=1:t
        lnl(i) = -n*0.5*log(2*pi) + log(abs(det(b))) - 0.5*log(det(omega)) - 0.5*u(i,:)*inv(omega)*u(i,:)';
    end
    
	rval = lnl;
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