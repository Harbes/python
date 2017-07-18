%=========================================================================
%
%   Likelihood ratio tests of deterministic trends in a vecm
%
%=========================================================================
function coint_trend( )

    clear all
    clc

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234) );

    t = 200; 
    n = 2;

    % Parameters: True DGP based on Model 3
    beta   = [1;-1]; 
    beta0  = 4; 
    beta1  = 0;
    alpha  = [-0.2;0.2]; 
    alpha0 = [1;1];
    delta0 = alpha0*0.2; 
    delta1 = alpha0*0;
    psi1   = [-.2 0;0 .4];

    % Simulate the model  
    v = randn(t,n);
    y = zeros(t+2,2);

    for j= 3:t+2 

        u = beta0 + beta1*j + y(j-1,:)*beta;        
        y(j,:) = y(j-1,:) + (delta0 + delta1*j + alpha*u)' ...
                 + (y(j-1,:) - y(j-2,:))*psi1 + v(j-2,:);
    end
    y   = trimr(y,2,0);
    mu0 = alpha*beta0 + delta0;
    mu1 = alpha*beta1 + delta1;

    % Estimate the vecm by maximum likelihood using Johansen estimator 
    p = 2;      % Number of lags in VAR      
    r = 1;      % Number of cointegrating equations    

    % Estimate Model 5    
    [~,~,~,lnl5,~,~] = johansen(y,p,r,5);  
    
    % Estimate Model 4    
    [~,~,~,lnl4,~,~] = johansen(y,p,r,4); 
    
    % Estimate Model 3     
    [~,~,~,lnl3,~,~] = johansen(y,p,r,3);     

    
    % LR test of Model 5 (unrestricted) against Model 4 (restricted)    
    lr = -2*t*(lnl4 - lnl5);

    disp(['Log-likelihood (unrestricted model 5) = ',num2str(lnl5) ]);
    disp(['Log-likelihood (restricted model 4)   = ',num2str(lnl4) ]);
    disp(['Value of likelihood ratio statistic   = ',num2str(lr) ]);
    disp(['p-value                               = ',num2str(1-chi2cdf(lr,n-r)) ]);
    disp(' ')

    % LR test of Model 4 (unrestricted) against Model 3 (restricted)    
    lr = -2*t*(lnl3 - lnl4);

    disp(['Log-likelihood (unrestricted model 4) = ',num2str(lnl4) ]);
    disp(['Log-likelihood (restricted model 3)   = ',num2str(lnl3) ]);
    disp(['Value of likelihood ratio statistic   = ',num2str(lr) ]);
    disp(['p-value                               = ',num2str(1-chi2cdf(lr,1)) ]);
    disp(' ')

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
    logl   = 0.5*(-tc*(log(2*pi)+1) - log(det(s00)) - sum( log(1-lam(1:r))));
    tracet = -tr*flipud(cumsum(log(1-flipud(lam))));                     
    maxt   = -tr*log(1-lam);                  
    
end

