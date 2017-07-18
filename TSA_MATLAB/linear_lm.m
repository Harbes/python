%=========================================================================
%
%    Program to perform a Lagrange multiplier test of a regression model
%
%=========================================================================
function linear_lr( ) 		

    clear all;
    clc;

    t = 200;
    %flag = 1;       % 1 = simulate data
    flag = 0;        % 0 = use GAUSS data to reproduce text
   
    if flag 
        [ y,x ] = simulatedata( t );
    else      
        % Aternatively load the GAUSS data
        load linear_testdata.dat
        y = linear_testdata(:,1);
        x = linear_testdata(:,[2 3 4]);     
    end
    x1 = x(:,2);
    x2 = x(:,3);
       
    % Estimate the constrained model  
    start = rand(3,1);
    theta = fminunc(@(theta) neglog0(theta,y,x1,x2),start);

% Lagrange Multiplier test (based on numerical opg matix) 
    theta0 = [theta([1 2]); (1-theta(2)); theta(3)] ;  	  
    gmat = numgrad(@lnl1,theta0,y,x1,x2);                               
    g    = mean(gmat)';
    j    = gmat'*gmat/t;
    lm   = t*g'*inv(j)*g;
    vcov = (1/t)*inv(j);


	disp('Gradient evaluated at contrained estimates');
    disp( g );

	disp('Outer product of gradients matrix');
    disp( j );
    
    disp('Covariance matrix');
    disp( vcov );

    disp(['LM test = ', num2str(lm) ]);
    disp(['p-value = ', num2str(1-cdf('chi2',lm,1) )]); 

end

%----------------------- Functions --------------------------%
% Simulate data
function [ y,x ] = simulatedata( t )

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123457) );

    beta0  = 1.0;                               % Parameter values                            
    beta1  = 0.7;
    beta2  = 0.3;
    sig   = sqrt(4.0);                              

    x1 = randn(t,1);                                                                 
    x2 = randn(t,1);
    u  = sig*randn(t,1);                         % Generate disturbances          
    y  = beta0 + beta1*x1 + beta2*x2 + u;                                           
    x  = [ones(t,1) x1 x2];

end 
% Constrained log-likelihood function
function lf = neglog0(theta,y,x1,x2)
    lf = -mean(lnl0(theta,y,x1,x2));
end
% Constrained log-likelihood function at each t    
function rval = lnl0(theta,y,x1,x2)
	m =  theta(1) + theta(2)*x1 + (1-theta(2))*x2;   % Mean                        
    s2 = theta(3);                                   % Variance                    
    z =  (y - m)/sqrt(s2);                           % Standardised residual       
    
    rval = -0.5*log(2*pi) - 0.5*log(s2) - 0.5*z.^2;
end

% Unconstrained log-likelihood function
function lf = neglog1(theta,y,x1,x2)
    lf = -mean(lnl1(theta,y,x1,x2));
end
% Unconstrained log-likelihood function at each t
        
function rval = lnl1(theta,y,x1,x2)
	m    = theta(1) + theta(2)*x1 + theta(3)*x2;       % Mean                        
    s2   = theta(4);                                   % Variance                    
    z    = (y - m)/sqrt(s2);                           % Standardised residual       
    
    rval = -0.5*log(2*pi) - 0.5*log(s2) - 0.5*z.^2;
    
end
