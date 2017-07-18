% ========================================================================
%
%      Program to estimate a simultaneous model with first order vector
%      autocorrelation. The set of equations is defined as yt*b + xt*a = u
%      where
%              u = ru(-1) + v
%          where yt is a (1xn) set of dependent variables at time t
%                xt is a (1xk) set of explanatory variables at time t
%
% ========================================================================
function auto_system( )

    clear all;
    clc;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123456) )
    
    t = 500;           
    
    % Simulate the data
    %[ y,x ] = simulatedata( t );
    
    % Load GAUSS data to reproduce results
    load system.mat
    
    y = gaussdata(:,[1 2]);
    x = gaussdata(:,[3 4]);
    
    % Estimate the unconstrained model  
    theta =  [0.6; 0.4; 0.2; -0.5; 0.8; 0.1;-0.2; 0.6];
    [theta1,a1,~,~,~,H] = fminunc(@(theta) neglog1(theta,y,x),theta);

    lnl1  = -(t-1)*a1;                           % Unconstrained log-likelihood                       
    vcov1 = inv(H);

    % Estimate the constrained model 
    theta =  [0.6; 0.4; 0.2; -0.5];
    [theta0,a0] = fminunc(@(theta) neglog0(theta,y,x),theta);

    lnl0 = -(t-1)*a0;                           % Constrained log-likelihood      

    % Estimate the constrained model with independent disturbances
    theta =  [0.6; 0.4; 0.2; -0.5; 0.8; 0.1];
    [theta2,a2] = fminunc(@(theta) neglog2(theta,y,x),theta);

    lnl2 = -(t-1)*a2;                            % Constrained log-likelihood 

    disp(' ');
    disp(['Unconstrained log-likelihood function                  = ', num2str(lnl1) ] );
    disp(['Constrained log-likelihood function                    = ', num2str(lnl0) ] );
    disp(['Constrained log-likelihood function (independent auto) = ', num2str(t*lnl2) ] );
	disp(' ');


    % LR test of no autocorrelation     
    lr  = -2*(lnl0 - lnl1);
    dof = length(theta1) - length(theta0);
    disp(['LR test (no autocorrelation)       = ', num2str(lr) ]);
    disp(['Degrees of freedom                 = ', num2str(dof) ]);
    disp(['p-value                            = ', num2str(1-chi2cdf(lr,dof)) ]);
	disp(' ');


    % LR  test of independent autocorrelation
    lr  = -2*(lnl2 - lnl1);
    dof = length(theta1) - length(theta2);
    disp(['LR test (indep autocorrelation)       = ', num2str(lr) ]);
    disp(['Degrees of freedom                    = ', num2str(dof) ]);
    disp(['p-value                               = ', num2str(1-chi2cdf(lr,dof)) ]);
	disp(' ');


    %     Wald test of no autocorrelation                       

    r  = [0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 ;
          0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 ;
          0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 ;
          0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 ];

    q  = [ 0 ; 0 ; 0 ; 0];
    
    wd = t*(r*theta1 - q)'*inv(r*vcov1*r')*(r*theta1 - q);
    dof = size(r,1);
    disp(' ');
    disp(['Wald test (no autocorrelation)       = ', num2str(wd) ]);
    disp(['Number of degrees of freedom         = ', num2str(dof) ]);
    disp(['p-value                              = ', num2str(1-chi2cdf(wd,dof)) ]);


    %     Wald test of common autocorrelation                   
    r  = [ 0 , 0 , 0 , 0 , 1 , 0 ,  0 , -1 ;
           0 , 0 , 0 , 0 , 0 , 1 , -1 ,  0 ];

    q  = [0 ; 0] ;
    wd = t*(r*theta1 - q)'*inv(r*vcov1*r')*(r*theta1 - q);
    dof = size(r,1);

    disp(' ');
    disp(['Wald test (indep autocorrelation)    = ', num2str(wd) ]);
    disp(['Number of degrees of freedom         = ', num2str(dof) ]);
    disp(['p-value                              = ', num2str(1-chi2cdf(wd,dof)) ]);


    % Lagrange Multiplier test (based on numerical opg matix)              
    dof    = length(theta1) - length(theta0);
    theta  = [theta0 ; zeros(4,1)];                  
    gmat   = numgrad(@lnlt1,theta,y,x);
    g      = mean(gmat)';
    j      = gmat'*gmat/t;
    lm     = t*g'*inv(j)*g;


    disp(' ');
    disp(['LM test (no autocorrelation)         = ', num2str(lm) ]);
    disp(['Number of degrees of freedom         = ', num2str(dof) ]);
    disp(['p-value                              = ', num2str(1-chi2cdf(lm,dof)) ]);

end
%
% -------------------------- Functions ---------------------------------
%
%-----------------------------------------------------------------------
%   Simulate the data  
%-----------------------------------------------------------------------
function [ y,x ] = simulatedata( t )

    %  Population paramaters     
    beta1  = 0.6; 
    alpha1 = 0.4;
    beta2  = 0.2; 
    alpha2 = -0.5; 

    rho11 = 0.8; 
    rho12 = 0.1;
    rho21 = -0.2; 
    rho22 = 0.6;

    omega = [ 1   0.5 ;
            0.5  1.0 ];
   
    b  =  [     1   -beta2 ; 
            -beta1     1   ]; 

    a  = [  -alpha1    0  ;
            0    -alpha2   ];

    % Exogenous variables                     
    x = [10*rand(t,1), 3*randn(t,1)];

    % Disturbances                          
    v = randn(t,2)*chol(omega);
    u = zeros(t,2);

    for i = 2:t;

        u(i,1) = rho11*u(i-1,1) + rho12*u(i-1,2) + v(i,1);
        u(i,2) = rho21*u(i-1,1) + rho22*u(i-1,2) + v(i,2);

    end
    
    % Simulate the model by simulating the reduced form
    y = zeros(t,2);

    i = 1;
    for i = 1:t
        y(i,:) = -x(i,:)*a*inv(b) + u(i,:)*inv(b);
    end   
end
%-----------------------------------------------------------------------
% Negative unconstrained log-likelihood  
%-----------------------------------------------------------------------
function lf = neglog1(theta,y,x)

    lf = -mean( lnlt1(theta,y,x) );

end
%-----------------------------------------------------------------------
% Unconstrained log-likelihood function
%-----------------------------------------------------------------------
function lnl = lnlt1(theta,y,x)

    [t,n] = size(y);
    b     = [    1    , -theta(3) ;
               -theta(1)  ,    1   ]; 
    a     = [ -theta(2) ,     0  ;
                 0     , -theta(4) ];

    rho   = [theta(5)  ,  theta(7) ;
             theta(6)  ,  theta(8) ];

    % Construct residuals and concentrate the covariance matrix
    u = zeros(t,n);
    v = zeros(t,n);
    for i = 2:t
        u(i,:) = y(i,:)*b + x(i,:)*a;     
        v(i,:) = u(i,:) - u(i-1,:)*rho;               
    end
    omega = v'*v/t;  
    
    lnl = zeros(t,1);
    for i = 2:t
        lnl(i) = - n*0.5*log(2*pi) + log(abs(det(b))) - 0.5*log(det(omega)) ...
            - 0.5*v(i,:)*inv(omega)*v(i,:)';
    end
    lnl =  trimr(lnl,1,0);                    
end
%-----------------------------------------------------------------------
% Negative constrained log-likelihood function
%-----------------------------------------------------------------------
function lf = neglog0(theta,y,x)

    lf = -mean( lnlt0(theta,y,x) );

end
%-----------------------------------------------------------------------
% Constrained log-likelihood function
%-----------------------------------------------------------------------  
function lnl = lnlt0(theta,y,x)

    [t,n] = size(y);
    b     = [    1    , -theta(3) ;
               -theta(1)  ,    1   ]; 
    a     = [ -theta(2) ,     0  ;
                 0     , -theta(4) ];
    rho   = [0  , 0 ;
             0  , 0 ];

    u = zeros(t,n);
    v = zeros(t,n);
    
    % Construct residuals and concentrate covariance matrix
    for i = 2:t

        u(i,:) = y(i,:)*b + x(i,:)*a;     
        v(i,:) = u(i,:) - u(i-1,:)*rho;                 
       
    end
    omega = v'*v/t;    
    
    lnl   = zeros(t,1);  
    for i = 2:t
        lnl(i) = - n*0.5*log(2*pi) + log(abs(det(b))) - 0.5*log(det(omega)) ...
            - 0.5*v(i,:)*inv(omega)*v(i,:)';
    end
    lnl = trimr(lnl,1,0);            
    
end
%-----------------------------------------------------------------------
% Negative constrained log-likelihood function (independent auto)
%-----------------------------------------------------------------------
function lf = neglog2(theta,y,x)

    lf = -mean( lnlt2(theta,y,x) );

end
%-----------------------------------------------------------------------
%   Constrained log-likelihood function at each observation
%   with independent autocorrelation
%-----------------------------------------------------------------------
function lnl = lnlt2(theta,y,x)

    [t,n] = size(y);
    b     = [     1    , -theta(3) ;
              -theta(1),    1      ]; 
    a     = [ -theta(2),     0     ;
                 0     , -theta(4) ];
    rho   = [ theta(5) ,     0     ; 
               0       ,  theta(6) ];

    % Construct residuals and concentrate covariance matrix
    u = zeros(t,n);
    v = zeros(t,n);
    for i = 2:t
        u(i,:) = y(i,:)*b + x(i,:)*a;     
        v(i,:) = u(i,:) - u(i-1,:)*rho;                    
    end
    omega = v'*v/t;
    
    lnl = zeros(t,1);
    for i = 2:t
        lnl(i) = - n*0.5*log(2*pi) + log(abs(det(b))) - 0.5*log(det(omega)) ...
            - 0.5*v(i,:)*inv(omega)*v(i,:)';
    end
    lnl = trimr(lnl,1,0);  
end


