%=========================================================================
%
%    Program to perform a Wald test of a regression model
%
%=========================================================================
function linear_wd( ) 		

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
    start = rand(4,1);                         
    theta = fminunc(@(theta) neglog1(theta,y,x1,x2),start); 
    ht    = numhess(@neglog1,theta,y,x1,x2);               % As f = negative log likelihood
    vcov1 = (1/t)*inv(ht);

    % Wald test       
    r = [0 1 1 0];
    q = 1;
    w = (r*theta - q)'*inv(r*vcov1*r')*(r*theta - q);


    fprintf('\nvcov based on numerical derivatives\n');
    for i = 1:size(vcov1,1)
        fprintf('%2.5f         %2.5f    %2.5f       %2.5f\n',vcov1(i,:));
    end 

    fprintf('\nWald test = %f\n', w);
    fprintf('p-value   = %f\n\n', 1-cdf('chi2',w,1) );

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
