%=========================================================================
%
%    Program to estimate linear regression by maximum likelihood.
%
%=========================================================================
function linear_estimate( ) 		

    clear all;
    clc;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) )

    t = 200;
    
    % Parameter values
    beta0  = 1.0;                                                         
    beta1  = 0.7;
    beta2  = 0.3;
    sig    = sqrt(4.0);
   
    % Simulate the data
    x1 = randn(t,1);                                                                 
    x2 = randn(t,1);
    u  = sig*randn(t,1);                                 
    y  = beta0 + beta1*x1 + beta2*x2 + u;                                       
    x  = [ones(t,1)  x1  x2];

    data = [y x1 x2 ];

    % Estimate the model
    theta_0 = rand(4,1);   % Initial guess 
    options = optimset('LargeScale','off','Display','iter');

    theta  = fminunc(@(theta) neglog(theta,y,x1,x2),theta_0,options);

    disp( 'Estimated Parameters')
    disp( theta )
    
    ht = numhess(@neglog,theta,y,x1,x2);
    disp( 'Estimated covariance matrix');
    disp((1/t)*inv(ht));
 

    % Estimate the concentrated model 
    theta  = fminunc(@(theta) neglogc(theta,y,x1,x2),rand(3,1));

    disp( 'Estimated Parameters (Concentrated)')
    disp( theta )
    
    htc = numhess(@neglogc,theta,y,x1,x2);
    disp( 'Estimated covariance matrix (concentrated)');
    disp((1/t)*inv(htc));
    

    %Compute OLS estimates
    theta = x\y;
    disp('OLS Parameter estimates')
    disp( theta );

    % Compute the covariance matrices    
    e       = y - x*theta;           
    sig2hat = e'*e/t;
    vcov    = sig2hat*inv(x'*x);
    
    disp('Covariance matrix (OLS)');
    disp( vcov );
    
    disp(['Variance estimate (OLS)          = ', num2str(sig2hat) ] );
    disp(['Variance estimate standard error = ', num2str(2*sig2hat^2/t) ]);


end

%
%--------------------------- Functions -----------------------------------
% 

%-------------------------------------------------------------------------
% Log-likelihood function 
%-------------------------------------------------------------------------

function lf = neglog(theta,y,x1,x2)
   
    lf = -mean( lnlt(theta,y,x1,x2) );

end

%-------------------------------------------------------------------------
% Log-likelihood function concentrated
%-------------------------------------------------------------------------

function lf = neglogc(theta,y,x1,x2)
  
    lf = -mean( lnltc(theta,y,x1,x2) );

end

%-------------------------------------------------------------------------
% Log-likelihood function at each observation
%-------------------------------------------------------------------------
             
function lf = lnlt(theta,y,x1,x2)
	
    m  = theta(1) + theta(2)*x1 + theta(3)*x2;   	                        
    s2 = theta(4);                              	                   
    z  = (y - m)/sqrt(s2);     
	lf = -0.5*log(2*pi) - 0.5*log(s2) - 0.5*z.^2;
end
    

%-------------------------------------------------------------------------
% Concentrated Log-likelihood function at each observation
%-------------------------------------------------------------------------
function lf = lnltc(theta,y,x1,x2)
	
	m  = theta(1) + theta(2)*x1 + theta(3)*x2;   	
    u  = y - m;                                  	                  
    s2 = u'*u/length(y);                                 		
	z  = (y - m)/sqrt(s2);    
	lf = -0.5*log(2*pi) - 0.5*log(s2) - 0.5*z.^2;
end



