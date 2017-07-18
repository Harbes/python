%=========================================================================
%
%    Program to perform a likelihood ratio test of a regression model
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
    
    % Estimate the unconstrained model  
    theta_0 = rand(4,1);                         
    [ theta1,lf1 ] = fminunc(@(theta) neglog1(theta,y,x1,x2),theta_0); 
    lf1 = -lf1;
    
    % Estimate the constrained model  
    theta_0 = rand(3,1);
    [ theta0,lf0 ] = fminunc(@(theta) neglog0(theta,y,x1,x2),theta_0);
    lf0 = - lf0;

    fprintf('\nParameter estimates (unconstrained)   = %10.6f %10.6f %10.6f %10.6f', theta1');
    fprintf('\nLog of the likelihood (unconstrained) = %f\n\n', lf1);
    fprintf('');

    %th0 = [theta(1:2,1) ; 1-theta(2); theta(3)];
    fprintf('Parameter estimates (constrained)   = %10.6f %10.6f %10.6f \n', theta0');
    fprintf('Log of the likelihood (constrained) = %f\n\n', lf0);

    % Likelihood ratio test       
    lr = -2*(t*lf0 - t*lf1);
    
    fprintf('Likelihood ratio test               = %f\n', lr);
    fprintf('p-value               			    = %f\n\n', 1-cdf('chi2',lr,1) );

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
