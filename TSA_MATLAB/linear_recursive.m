% ========================================================================
%
%   Program to estimate a recursive system to demonstrate the 
%   equivalence of OLS and MLE in this special case of the SUR model.
%
%   The model is  
%       y1t =                       + alpha1*x1t + u1t
%       y2t = beta1*y1t             + alpha2*x2t + u2t
%       y3t = beta2*y1t + beta3*y2t + alpha3*x3t + u3t
%
% ========================================================================

function linear_recursive(  )

    clear all;
    clc;

    % Simulate data
    t = 200;   
    [ y,x ] = simulatedata( t );
    
    % Estimate the model with random starting values
    theta = fminunc( @(theta) neglog(theta,y,x),rand(6,1) );
    
    % Estimate parameters by OLS
    eq1 = x(:,1)\y(:,1);       
	eq2 = [y(:,1) x(:,1)]\ y(:,2);
	eq3 = [y(:,[1 2]) x(:,1)]\y(:,3);
  
    true     = [0.4; 0.6; -0.5; 0.2; 1.0; 0.2];
    thetaols = [eq1; eq2; eq3 ];

    disp('Comparing true and estimated parameter values')
    disp('    Actual    MLE     OLS')
    disp( [ true theta thetaols] );

    % Compute covariance matrix
    u = zeros(t,3);
    b = [1 -theta(2) -theta(4); 0  1  -theta(5); 0  0 1];
    a = [-theta(1) -theta(3) -theta(6); 0  0  0; 0  0  0];

    for i=1:t
        u(i,:) = y(i,:)*b + x(i,:)*a;     	                    
    end
    tmp    = diag(u'*u/t);
    omegah = eye(3);
    omegah = diag( tmp );
    omega  = [2 0 0; 0 1 0; 0 0 5];

    disp('Comparison of true and estimated elements of omega');
    disp( [ omega(:) omegah(:) ] );

end

%
%--------------------- Functions ----------------------%
%
% Simulate data
function [ y,x ] = simulatedata( t )

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123457) );

    alpha1 = 0.4;
    beta1  = 0.6;
    beta2  = 0.2; 
    alpha2 = -0.5;  
    beta3  = 1.0; 
    alpha3 =  0.2; 
    omega  = [2 0 0; 0 1 0; 0 0 5];


    % Construct population parameter matrices     
    b = [1 -beta1 -beta2; 0  1  -beta3; 0  0  1];
    a = [-alpha1 -alpha2 -alpha3; 0  0  0; 0  0  0];

    % Construct exogenous variables                       
    x = [randn(t,1) 2*randn(t,1) 3*randn(t,1)];

    % Construct disturbances                              
    u = randn(t,3)*chol(omega);

    % Simulate the model by simulating the reduced form   
    y = zeros(t,3);

    for i=1:t
        y(i,:) = -x(i,:)*a*inv(b) + u(i,:)*inv(b);
    end
end

% Negative log-likelihood function
function lf = neglog( theta,y,x )

    lf = -mean( lnlt(theta,y,x) );

end

% Log-likelihood at each observation         
function lnl = lnlt(theta,y,x,u)
	
    [t n] = size(y); 
    b = [1 -theta(2) -theta(4); 0  1  -theta(5); 0  0 1];
    a = [-theta(1) -theta(3) -theta(6); 0  0  0; 0  0  0];
            
    % Construct residuals and concentrate diagonal covariance matrix
    u = zeros(t,n);
    for i=1:t
        u(i,:) = y(i,:)*b + x(i,:)*a;     	                    
    end
    tmp   = diag(u'*u/t);
    omega = eye(n);
    omega = diag( tmp );

    lnl = zeros(t,1);   
   	for i=1:t
	    lnl(i) = -n*0.5*log(2*pi) + log(abs(det(b))) - 0.5*log(det(omega)) ...
            - 0.5*u(i,:)*inv(omega)*u(i,:)';
    end

end
	
