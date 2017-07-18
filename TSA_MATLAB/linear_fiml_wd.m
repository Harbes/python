%=========================================================================
%
%  Program to estimate model by full information maximum likelihood 
%  and do a Wald test. 
%
%=========================================================================
function linear_fiml_wd( ) 		 

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
    lf1   = -lf1;
    ht    = numhess(@neglog1,theta1,y,x);               
    vcov1 = (1/t)*inv(ht);

    fprintf('\n\nParameter estimates (unconstrained)   = \n%10.6f %10.6f\n%10.6f %10.6f\n', theta1');
    fprintf('Log of the likelihood (unconstrained) = %f\n\n', lf1);

    % Wald test
    q = 0;
    r = [0 1 0 1];
    w = (r*theta1 - q)'*inv(r*vcov1*r')*(r*theta1 - q);


    fprintf('\nCovariance matrix\n');
    for i = 1:size(vcov1,1)
        fprintf('%10.6f %10.6f %10.6f %10.6f\n',vcov1(i,:)*1000);
    end 

    fprintf('\nWald test   = %f\n', w);
    fprintf('p-value     = %f\n\n', 1-cdf('chi2',w,1) );


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
