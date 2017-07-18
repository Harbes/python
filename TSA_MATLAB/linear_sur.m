%=========================================================================
%
%     Program to estimate a SUR system
%
%=========================================================================

function linear_sur( ) 		


    clear all;
    clc;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) )


    t = 500;                

    alpha1 =  0.4;
    alpha2 = -0.5;
    alpha3 =  1.0;

    omega = [   1		0.5		-0.1;
                0.5 	1.0 	 0.2;
                -0.1 	0.2 	 1.0];    

    b = eye(3); 
    a = [	-alpha1 	0       	0;
        0    		-alpha2     0;
        0        	0      		-alpha3];


    % Exogenous variables                       
    x = [randn(t,1) 2*randn(t,1) 3*randn(t,1)];


    % Disturbances                              
    u = randn(t,3)*chol(omega);


    % Simulate the model by simulating the reduced form   
    % Note that the SUR system is the reduced form

    y = zeros(t,3);

    for i=1:t
        
    y(i,:) = -x(i,:)*a*inv(b) + u(i,:)*inv(b);
    
    end   

    % Estimate the model  
    theta0 = rand(3,1);
    theta  = fminunc(@(theta) neglog(theta,y,x),theta0);

    % Compute OLS estimates for comparative purposes
    alpha1ols =  x(:,1)\y(:,1);
    alpha2ols =  x(:,2)\y(:,2);
    alpha3ols =  x(:,3)\y(:,3);
    
    alpha    = [alpha1; alpha2; alpha3];
    alphaols = [alpha1ols; alpha2ols; alpha3ols];
    
    disp('Comparing true and estimated parameter values')
    disp('    Actual    MLE       OLS')
    disp( [ alpha theta alphaols] );

    % Compute residuals at optimal parameters
    a     = [	-theta(1) 	0        	0;
                0      		-theta(2)   0;
                0        	0        	-theta(3) ];          

    for i=1:t;
        
        u(i,:) = y(i,:)*b + x(i,:)*a ; 
        
    end
    omegahat = u'*u/t;

    disp('Comparing true and estimated elements of Omega') 
    disp( [ omega(:) omegahat(:) ] );

end

%
%------------------------- Functions -------------------------------------%
%

% Log-likelihood function

function lf = neglog( theta,y,x )

    lf = -mean( lnlt( theta,y,x ) );

end

% Log-likelihood function at each observation     

function lnl = lnlt( theta,y,x ) 
	
    
    [t n] = size(y);  
    b     = eye(n); 
    a     = [	-theta(1) 	0        	0;
                0      		-theta(2)   0;
                0        	0        	-theta(3) ];          
    u     = zeros(t,n);
    
    % Compute residuals	
    for i=1:t;
        
        u(i,:) = y(i,:)*b + x(i,:)*a; 
        
    end
   
    % Concentrate residual var-covar matrix  
    omega = u'*u/t;                       	
    
    % Log-likelihood function
    lnl = zeros(t,1);
    
   	for i=1:t
        
	    lnl(i) = -n*0.5*log(2*pi) + log(abs(det(b))) - 0.5*log(det(omega)) - 0.5*u(i,:)*inv(omega)*u(i,:)';
        
    end
    
	
end


