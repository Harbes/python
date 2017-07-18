%=========================================================================
%
%   Program to estimate bivariate term structure model
%
%=========================================================================
function coint_bivterm( )

    clear all
    clc

    % Load data and choose 10-year and 1-year yields
    load usmacro.mat   
    yield = [r10yr  r1yr];
    t     = length(yield);
    
    z0   = trimr(yield,1,0) - trimr(yield,0,1);
    z1   = trimr(yield,0,1);
    
    % Full rank case: Estimate the var in levels 
    x = [ ones(t-1,1)  trimr(yield,0,1) ];
    y = trimr(yield,1,0);
    [k,n]  = size(y);
    b      = x\y;
    v      = y - x*b;
    omegav = v'*v/t;

    lnl = -0.5*k*n*(log(2*pi) + 1) - 0.5*k*log(det(omegav));
    disp('Estimated coefficients');
    disp( b )
    disp(' ' );
    disp('Estimated covariance matrix');
    disp( omegav );

    disp(['Log-likelihood of full rank model    = ', num2str(lnl/t) ]);
    disp(['Determinant of omegav                = ', num2str(det(omegav))]);

    % Reduced rank case parameters
    p     = 1;      
    r     = 1;    
    model = 2;

    % Reduced rank case: iterative procedure (Model 2)
    % Starting values
    b = [ones(t,1) yield(:,2)]\yield(:,1);
    u = yield(:,1) - [ ones(t,1) yield(:,2)]*b;
    a = trimr(u,0,1)\z0;

    pstart = [b ; a'];
    ops    = optimset('LargeScale', 'off', 'Display', 'iter');     
    [theta,logl,aa,aa,aa,hess] = fminunc(@(p) neglog(p,z0,z1),pstart,ops);

    disp('Iterative procedure estimates')        
    disp('-----------------------------')
    disp(['Log likelihood function = ',num2str(-logl) ])
    disp(['beta_c         = ',num2str(theta(1)) ])
    disp(['beta_r         = ',num2str(theta(2)) ])
    disp(['alpha_1        = ',num2str(theta(3)) ])
    disp(['alpha_2        = ',num2str(theta(4)) ])
    vc = (1/t)*inv(hess);
    disp('Covariance matrix (iterative estimator)')
    disp(vc)
    disp(' ')
    
    % Perform Wald test of (1,-1) cointegrating vector)        
    wd = (theta(2)-1)^2/vc(2,2);      
    
    disp('Wald test of (1,-1) cointegrating vector');
    disp(['Wald statistic           = ',num2str(wd) ]);
    disp(['p-value                  = ',num2str(1-chi2cdf(wd,1)) ]);
    disp(' ')

    % Perform Wald test of y1 weakly exogenous    
    wd = (theta(3)-0)^2/vc(3,3);                         

    disp('Wald test of y1 weakly exogenous');
    disp(['Wald statistic           = ',num2str(wd) ]);
    disp(['p-value                  = ',num2str(1-chi2cdf(wd,1)) ]);
    disp(' ')
    
    % Perform Wald test of y2 weakly exogenous     
    wd = (theta(4)-0)^2/vc(4,4); 
    
    disp('Wald test of y2 weakly exogenous');
    disp(['Wald statistic           = ',num2str(wd) ]);
    disp(['p-value                  = ',num2str(1-chi2cdf(wd,1)) ]);
    disp(' ')

    % Johansen estimator of reduced rank model
    [ alpha,beta,logl,maxt,tracet ] = johansen(y,p,r,model);

    disp('Johansen procedure estimates')
    disp('----------------------------')
    disp(['Log likelihood function = ',num2str(logl) ])
    disp(['beta_c         = ',num2str(-beta(3)) ])
    disp(['beta_r         = ',num2str(-beta(2)) ])
    disp(['alpha_1        = ',num2str(alpha(1)) ])
    disp(['alpha_2        = ',num2str(alpha(2)) ])
    

    % Zero rank case: VAR in first differences         
    x  = ones(length(yield)-1,1);
    y  = trimr(yield,1,0) - trimr(yield,0,1);
    b  = x\y;
    v  = y;             
    vc = v'*v/length(v);
    n  = size(y,2);
    lf = -0.5*n*(log(2*pi) + 1) - 0.5*log(det(vc));

    disp(['Log-likelihood of zero rank model    = ',num2str(lf) ]);
    disp(['Determinant of covariance matrix     = ',num2str(det(vc))]);
    disp('Covariance matrix of residuals') 
    disp(vc)



    
end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Johansen procedure 
%-------------------------------------------------------------------------
function [alpha,beta,logl,maxt,tracet] = johansen(y,p,r,model)

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
    
    % Statistics
    logl   = -0.5*tc*(log(2*pi)+1) - 0.5*log(det(s00)) - 0.5*sum( log(1-lam(1:r)));
    tracet = -tr*flipud(cumsum(log(1-flipud(lam))));                     
    maxt   = -tr*log(1-lam);                  
    
end
%-------------------------------------------------------------------------
%  Reduced rank log-likelihood function: Model 2 and p=1
%-------------------------------------------------------------------------
function lf = neglog(b,z0,z1)

    nobs = length(z0);
    f   = zeros(nobs,1);

    v1 = z0(:,1) - b(3)*(z1(:,1) - b(1) - b(2)*z1(:,2));
    v2 = z0(:,2) - b(4)*(z1(:,1) - b(1) - b(2)*z1(:,2));
    v  = [v1  v2];

    [k n] = size( v );
    
    omegav = v'*v/k;

    for t = 1:nobs

        f(t) = - n*0.5*log(2*pi) - 0.5*log(det(omegav)) - 0.5*v(t,:)*inv(omegav)*v(t,:)';

    end

    lf = - mean( f );
end


