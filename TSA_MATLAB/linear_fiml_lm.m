%=========================================================================
%
%  Program to estimate model by full information maximum likelihood 
%  and do a LM test. 
%
%=========================================================================
function linear_fiml_lm( ) 		 

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

    % Estimate the constrained model  
    start = rand(3,1);
    theta = fminunc(@(theta) neglog0(theta,y,x),start);

% Lagrange Multiplier test (based on numerical opg matix) 
    theta0 = [theta([1 2 3]); -theta(2)] ;  	  
    gmat = numgrad(@lnlt1,theta0,y,x);                               
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
