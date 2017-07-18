%=========================================================================
%
%   Program to estimate a bivariate vecm of the permanent income hypothesis
%
%=========================================================================
function coint_permincome( )

    clear all
    clc

    % Load data 
    load permincome.mat
    
   
    % Select desired sample (1984Q1 to 2005Q4)
    rc = rcpc(149:236);
    ry = rypc(149:236);
    y  = [ log(rc)  log(ry) ];
    
    % Reduced rank case: Johansen estimator (Model 3)      
    p     = 4;     
    r     = 1;     
    model = 3;  

    [ alpha,beta,param,logl,maxt,tracet ] = johansen(y,p,r,model);

    disp(['Estimates of cointegrating vector   = ', num2str(beta') ]);
    
    mu0 = param(r+1,:)';
    disp(['Estimates of alpha                  = ', num2str(alpha') ]);
    disp(['Estimates of mu0                    = ', num2str(mu0') ]);
    
    % Compute orthonormal complement of alpha
    alpha0 = null(alpha');      
    disp(['Orthogonal complement matrix alpha0 = ', num2str(alpha0') ]);

    % Computer constant term in cointegrating regression
    beta0 = (alpha'*alpha)\(alpha'*mu0);
    
    disp(['Estimate of beta0                   = ', num2str(beta0) ]);

  
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

    % Estimates
    beta  = rref(gam(:,1:r)')';                                               
    alpha = ((r1*beta)\r0)';    
   
    % Estimates of VECM
    tmpx  = [ z1*beta z2];
    param = tmpx\z0;
    
    % Statistics
    logl   = 0.5*(-tr*tc*(log(2*pi)+1) - tr*log(det(s00)) - tr*sum( log(1-lam(1:r))));
    tracet = -tr*flipud(cumsum(log(1-flipud(lam))));                     
    maxt   = -tr*log(1-lam);                  
    
end






