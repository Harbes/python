%=========================================================================
%
%     Program to estimate a simultaneous model with vector heteroskedasticity 
%     The system is defined as yt*b + xt*a = u      
%     where
%       u ~ N(0,V)
%    and
%       yt is a (1xn) set of dependent variables at time t
%       xt is a (1xm) set of mean explanatory variables at time t
%       wt is a (1xs) set of variance explanatory variables at time t
%
%=========================================================================
function hetero_system( )

    clear all;
    clc;


    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123456) )

    % Simulate the model   
    t = 2000;
    [y,x,w] = simulatedata( t );

    % Estimate the model
    theta0 = [ 0.6; 0.4; 0.2;-0.5; 1.0; 0.5; 0.5; 0.2; 2.0; 0.2 ];
    [theta,a,~,~,~,H] = fminunc(@(theta) neglog(theta,y,x,w),theta0);

    vcov = inv(H);
    disp(['Log-likelihood function     = ',num2str(-a) ]);
    disp( 'Parameter estimates and standard errors' );
    disp([ theta sqrt(diag(vcov)) ]);

    % Wald test of no vector heteroskedasticty                        
    r   = [0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 ;
           0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 ;
           0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 ];
    q   = [0 ; 0 ; 0 ];  
    wd  = t*(r*theta - q)'*inv(r*vcov*r')*(r*theta - q);
    dof = size(r,1);
    disp(['Wald statistic          = ',num2str(wd) ]); 
    disp(['p-value                 = ',num2str(1-cdf('chi2',wd,dof)) ]);


end
%
% -------------------------- Functions ---------------------------------
%
%-----------------------------------------------------------------------
%   Simulate the data  
%-----------------------------------------------------------------------
function [ y,x,w ] = simulatedata( t )

    beta1  = 0.6; 
    alpha1 = 0.4;
    beta2  = 0.2; 
    alpha2 = -0.5; 

    c11 = 1.0;
    c21 = 0.5;
    c22 = 2.0;

    d11  = 0.5;
    d21  = 0.2;
    d22  = 0.2;

    b  =  [  1   , -beta2 ; 
            -beta1 ,    1   ]; 
    a  = [-alpha1 ,   0     ;
            0    , -alpha2 ];
    c  =  [c11    ,   0     ;
            c21    ,  c22    ];
    d  =  [d11    ,   0     ;
            d21    ,  d22   ] ;

    % Exogenous variables                      
    x = [10*rand(t,1) , 3*randn(t,1)];         	   
    w = rand(t,1);                            

    % Disturbances                              
    u = zeros(t,2);
    for i = 1:t
        l      = c + d*w(i);
        u(i,:) = randn(1,2)*l';
    end

    % Simulate the reduced form   
    y = zeros(t,2);

    for i = 1:t
        y(i,:) = -x(i,:)*a*inv(b) + u(i,:)*inv(b);

    end
end
%-----------------------------------------------------------------------
% Negative unconstrained log-likelihood  
%-----------------------------------------------------------------------
function lf = neglog(theta,y,x,w)

    lf = -mean( lnlt(theta,y,x,w) );

end
%-----------------------------------------------------------------------
% Unconstrained log-likelihood function at each observation
%-----------------------------------------------------------------------
function lnl = lnlt(theta,y,x,w)

    t = size(y,1);
    n = size(y,2);
    b  =   [   1     , -theta(3) ;
          -theta(1) ,     1     ];     
    a  = [-theta(2) ,     0     ;
              0     , -theta(4)] ;
    c  =   [theta(5) ,     0     ;
           theta(7) ,  theta(9) ];
    d  =  [theta(6) ,     0     ;
           theta(8) ,  theta(10)];
       
    u   = zeros(t,n);
    lnl = zeros(t,1);
    
    for i=1:t
        
        u(i,:) = y(i,:)*b + x(i,:)*a;     
        l      = c + d*w(i);
        V      = l*l';
        lnl(i) = - n*0.5*log(2*pi) + log(abs(det(b))) - 0.5*log(det(V)) ...
            - 0.5*u(i,:)*inv(V)*u(i,:)';
     
     end

end
