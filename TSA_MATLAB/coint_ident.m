%=========================================================================
%
%   Program to investigate alternative identification strategies
%
%=========================================================================
function coint_ident( )

    clear all
    clc

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) ); 

    t = 200; 
    n = 3;

    % Parameters of true DGP based on Model 3   
    beta   = [  1  0  ; 0   1  ; -0.8   -0.8  ];          
    beta0  = [4 1]; 
    beta1  = [0  0];
    alpha  = [-0.2  0 ; 0  -0.3;  0.1      0.1]; 
    alpha0 = null(alpha');
    delta0 = alpha0*0.1; 
    delta1 = alpha0*0;
    psi1   = [-0.2 0  0; 0  0.4 0 ; 0 0  0.1]; 


    % Simulate the model: DGP is a vecm based on Model 3  
    v = randn(t,n);
    y = zeros(t+2,n);

   for j= 3:t+2 

        u = beta0 + beta1*j + y(j-1,:)*beta;        
        y(j,:) = y(j-1,:) + (delta0 + delta1*j + alpha*u')' ...
                 + (y(j-1,:) - y(j-2,:))*psi1 + v(j-2,:);
    end
    y   = trimr(y,2,0);
    %mu0 = alpha*beta0 + delta0;
    %mu1 = alpha*beta1 + delta1;


    % Estimate the vecm by maximum likelihood using Johansen estimator 
    p = 2;      % Number of lags in VAR      
    r = 1;      % Number of cointegrating equations    
    model = 3;   

    [alpha,beta,param,logl,maxt,tracet] = johansen(y,p,r,model);


    % Estimate the vecm by maximum likelihood using the iterative estimator for triangular cross-equation normalizations  **/
    dy = trimr(y,1,0)-trimr(y,0,1); 
    z0 = trimr(dy,p-1,0);
    z1 = trimr(y,p-1,1);

    z2 = [];
    
    for j = 1:p-1
        z2 = [ z2 trimr(dy,p-1-j,j) ];
    end
    
    % Model 3    
    z2 = [ones(length(y)-p,1)   z2];   
    nobs = length(z0);

    theta_0 = 0.1*ones(19,1);
    ops    = optimset('LargeScale', 'off', 'Display', 'off');     
    [theta,logl,~,~,~,hess] = fminunc(@(p) neglog(p,z0,z1,z2),theta_0,ops);
    
    disp('Cointegrating vectors (triangularization)')
    disp('with cross-equation restrictions')
    disp( [eye(2) ; [-theta(10) -theta(10) ]] );


end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Johansen procedure 
%-------------------------------------------------------------------------
function [alpha,beta,param,logl,maxt,tracet] = johansen(y,p,r,model)

    ty = length(y);
    dy = trimr(y,1,0)-trimr(y,0,1); 
    z0 = trimr(dy,p-1,0);
    z1 = trimr(y,p-1,1);
                  
    z2 = [];
    
    for j =1:p-1
         
        z2 = [ z2 trimr(dy,p-1-j,j) ];
    end
    
    if model == 1;
                                      
        z1 = z1;
        z2 = z2;
        
    elseif model == 2 
            
        z1 = [ trimr(y,p-1,1)  ones(ty-p,1) ];
        z2 = z2;
        
    elseif model == 3 
            
        z1 = z1;
        z2 = [ ones(ty-p,1)  z2 ];   
     
    elseif model == 4
            
        z1 = [z1  seqa(1,1,ty-p)'];
        z2 = [ ones(ty-p,1)  z2 ];   

    elseif model == 5
            
        z1 = z1;
        z2 = [ones(ty-p,1)  seqa(1,1,ty-p)'  z2];   
    end
        
    if p == 1 && model <= 2
        
         r0 = z0; 
         r1 = z1;
         
    else
        
        r0 = z0-z2*(z2\z0); 
        r1 = z1-z2*(z2\z1);
    end
          
    [ tr,tc ]  = size( r0 );
    
    % Construct sums of squares matrices
    s00 = r0'*r0/tr;                                        
    s11 = r1'*r1/tr;
    s01 = r0'*r1/tr;
    s10 = s01';
    
    % Solve eigenvalue problem
    l       = chol(s11)';                                               
    [e,tmp] = eig( inv(l)*s10*inv(s00)*s01*inv(l') );
       
   
    % Sort eigenvalues and store index
    [ lam,IX ] = sort( diag(tmp),'descend' ); 

    % Sort eigenvectors on eigenvalue index and normalize
    gam = (inv(l')*e);
    gam = gam(:,IX);
    
    tmp = gam(:,1:r);
    disp('Cointegrating vectors (nonnormalizaed)')
    disp(tmp);
    
    disp('Cointegrating vectors (diagonalized)')
    disp(tmp./diag(tmp(1:r,1:r)));

    % Estimates
    beta  = rref(gam(:,1:r)')';    
    
    disp('Cointegrating vectors (triangularization)')
    disp(beta);

    alpha = ((r1*beta)\r0)';    
   
    % Estimates of VECM
    tmpx  = [ z1*beta z2];
    param = tmpx\z0;
    
    % Statistics
    logl   = 0.5*(-tc*(log(2*pi)+1) - log(det(s00)) - sum( log(1-lam(1:r))));
    tracet = -tr*flipud(cumsum(log(1-flipud(lam))));                     
    maxt   = -tr*log(1-lam);                  
    
end
%-------------------------------------------------------------------------
%  Reduced rank log-likelihood function: Model 2 and p=1
%-------------------------------------------------------------------------
function lf = neglog(b,z0,z1,z2)

    nobs = length(z0);
    f   = zeros(nobs,1);

     m1  = b(1) + b(4)*(z1(:,1) - b(10)*z1(:,3)) + b(7)*(z1(:,2) - b(10)*z1(:,3)) + z2(:,[2 3 4])*b([11 12 13]);
     m2  = b(2) + b(5)*(z1(:,1) - b(10)*z1(:,3)) + b(8)*(z1(:,2) - b(10)*z1(:,3)) + z2(:,[2 3 4])*b([14 15 16]);
     m3  = b(3) + b(6)*(z1(:,1) - b(10)*z1(:,3)) + b(9)*(z1(:,2) - b(10)*z1(:,3)) + z2(:,[2 3 4])*b([17 18 19]);

     v1 = z0(:,1) - m1;
     v2 = z0(:,2) - m2;
     v3 = z0(:,3) - m3;

     v = [v1  v2  v3 ];

    [k n] = size( v );
    
    omegav = v'*v/k;

    for t = 1:nobs

        f(t) = - n*0.5*log(2*pi) - 0.5*log(det(omegav)) - 0.5*v(t,:)*inv(omegav)*v(t,:)';

    end

    lf = - mean( f );
end
