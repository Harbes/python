%=========================================================================
%
%  Program to estimate model by full information maximum likelihood 
%  and do a LR test. 
%
%=========================================================================
function linear_fiml_lr( ) 		 

    clear all;
    clc;
    t = 500;
    %flag = 1;       % 1 = simulate data
    flag = 0;        % 0 = use GAUSS data to reproduce text
   
    if flag 
        [ y,x ] = simulatedata( t );
    else      
        % Aternatively load the GAUSS data
        load linear_fimltestdata.dat
        y = linear_fimltestdata(:,[1 2]);
        x = linear_fimltestdata(:,[3 4]);     
    end

    % Estimate the unconstrained model  
    start = rand(4,1);                         
    [ theta1,lf1 ] = fminunc(@(theta) neglog1(theta,y,x),start); 
    lf1 = -lf1;
    u       = zeros(t,2);
    for i=1:t
        u(i,:) = y(i,:)*[1 -theta1(3); -theta1(1) 1] + x(i,:)*[-theta1(2) 0; 0 -theta1(4)];     	                    
    end

    disp('Residual variance-covariance matrix (unrestricted)')
    omega1 = u'*u/t;
    disp(omega1);

    disp(u'*u/t);
    
    % Estimate the constrained model  
    start = rand(3,1);
    [ theta0,lf0 ] = fminunc(@(theta) neglog0(theta,y,x),start);
    lf0 = - lf0;

    u       = zeros(t,2);
    for i=1:t
        u(i,:) = y(i,:)*[1 -theta0(3); -theta0(1) 1] + x(i,:)*[-theta0(2) 0; 0  theta0(2)];     	                    
    end
    disp('Residual variance-covariance matrix (restricted)')
    omega0 = u'*u/t;
    disp(omega0);

    
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

    % Alternative form of the LR test
    b1 = [1 -theta1(3); -theta1(1) 1];
    b0 = [1 -theta0(3); -theta0(1) 1];
    
    lr_alternative = t*( log(det(omega0)) - log(det(omega1)) ) - 2*t*( log(abs(det(b0))) - log(abs(det(b1))) );


    fprintf('Residual variance-covariance (unconstrained) = \n%10.6f %10.6f\n%10.6f %10.6f\n', omega1);
    fprintf('Residual variance-covariance (constrained)   = \n%10.6f %10.6f\n%10.6f %10.6f\n', omega0);

    fprintf('Likelihood ratio test (alternative)          = %f\n', lr_alternative);
    fprintf('p-value                                      = %f\n\n', 1-cdf('chi2',lr_alternative,1) );


end

%--------------------- Functions ----------------------%
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

% Unconstrained log-likelihood function            
function lf = neglog1(theta,y,x)
    
    lf = -mean(lnlt1(theta,y,x));
    
end

% Unconstrained log-likelihood function at each observation            
function lf = lnlt1(theta,y,x)	
    
    [t n] = size(y);
    
    b = [1 -theta(3); -theta(1) 1]; 
    a = [-theta(2) 0; 0 -theta(4)];
    u = zeros(t,n);

    for i=1:t
        u(i,:) = y(i,:)*b + x(i,:)*a;     	                    
    end

    omega = u'*u/t;                         
    lf    = zeros(t,1);

    for i=1:t
        lf(i) = -n*0.5*log(2*pi) + log(abs(det(b))) - 0.5*log(det(omega)) - 0.5*u(i,:)*inv(omega)*u(i,:)';
    end

end

% Constrained log-likelihood function            
function lf = neglog0(theta,y,x)
    
    lf = -mean(lnlt0(theta,y,x));
    
end


% Constrained log-likelihood function at each observation            
function lf = lnlt0(theta,y,x)	
    
    [t n] = size(y);
    
    b = [1 -theta(3); -theta(1) 1]; 
    a = [-theta(2) 0; 0 theta(2)];
    u = zeros(t,n);

    for i=1:t
        u(i,:) = y(i,:)*b + x(i,:)*a;     	                    
    end

    omega = u'*u/t;                         
    lf    = zeros(t,1);

    for i=1:t
        lf(i) = -n*0.5*log(2*pi) + log(abs(det(b))) - 0.5*log(det(omega)) - 0.5*u(i,:)*inv(omega)*u(i,:)';
    end

end

