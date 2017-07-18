%=========================================================================
%
%   Program to test for weak exogeneity using Wald and LR tests
%
%=========================================================================
function coint_exogeneity( )

    clear all
    clc
    
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',1) );

    t = 200; 
    n=2;

    % Parameters: true DGP based on Model 3   
    beta      = [1; -1]; 
    beta0     = 0; 
    beta1     = 0;
    alpha     = [-0.4 ; 0]; 
    alphaorth = [ 0 ; 1 ];
    delta0    = alphaorth*0; 
    delta1    = alphaorth*0;
    psi1      = [-0.2 0 ; 0 0.4];
    omegav    = [1 0.5; 0.5 1];

    % Simulate the model   
    v = randn(t,2)*chol(omegav);
    y = zeros(t+2,2);

    for j = 3:t+2 

        u      = beta0 + beta1*j+y(j-1,:)*beta;
        y(j,:) = y(j-1,:) + (delta0 + delta1*j + alpha*u')' + (y(j-1,:) - y(j-2,:))*psi1 + v(j-2,:);

    end
    y = trimr(y,2,0);

    % Estimate the vecm by maximum likelihood  
    p     = 2;      % Number of lags in VAR      
    r     = 1;      % Number of cointegrating equations   
    model = 1;
    Tp = t-p;
    
    % Johansen estimator of reduced rank model
    [ ~,beta,~,~,~ ] = johansen(y,p,r,model);

    disp('Johansen estimate of long-run parameters')
    disp(beta)
    disp( ' ')

    
    dy = trimr(y,1,0)-trimr(y,0,1); 
    z0 = trimr(dy,p-1,0);
    z1 = trimr(y,p-1,1);

    z2 = [];
    for j = 1:p-1 
        
        z2 = [ z2 trimr(dy,p-1-j,j) ];
    end

    nobs = length(z0);

    % Estimate the unconstrained model
    theta_0                 = 0.1*ones(7,1);
    ops                     = optimset('LargeScale', 'off', 'Display', 'off');     
    [theta1,lf1,~,~,~,hess] = fminunc(@(p) neglog1(p,z0,z1,z2),theta_0,ops);

    lf1 = -lf1;
    vc  = (1/nobs)*inv(hess);

    % Wald test (y2 weakly exogenous)     
    r = [ 0   1   0   0   0   0   0 ];
    q = 0;
    w    = (r*theta1 - q)'*inv(r*vc*r')*(r*theta1 - q);

    disp('Wald test (y2 weakly exogenous)');
    disp(['Wald statistic           = ',num2str(w) ]);
    disp(['p-value                  = ',num2str(1-chi2cdf(w,1)) ]);
    disp(' ')


    % Wald test (y2 strongly exogenous)       
    r = [0   1   0   0   0   1   0];
    q = 0;
    w = (r*theta1 - q)'*inv(r*vc*r')*(r*theta1 - q);
    
    disp('Wald test (y2 strongly exogenous)');
    disp(['Wald statistic           = ',num2str(w) ]);
    disp(['p-value                  = ',num2str(1-chi2cdf(w,2)) ]);
    disp(' ')


    % Estimate constrained model (y2 is weakly exogenous)   
    theta_0             = [theta1(1) ; theta1(3:7)];
    [~,lf0,~,~,~,~] = fminunc(@(p) neglog0(p,z0,z1,z2),theta_0,ops);

    lf0 = -lf0;

    % LR test of weak exogeneity
    lr = -2*nobs*(lf0 - lf1);
    
    disp('LR test (y2 weakly exogenous)');
    disp(['LR statistic           = ',num2str(lr) ]);
    disp(['p-value                = ',num2str(1-chi2cdf(lr,1)) ]);
    disp(' ')

    % Estimate constrained model (y2 is strongly exogenous)   
    theta_0         = [theta1(1) ; theta1(3:5); theta1(7)];
    [~,lf00,~,~,~,~] = fminunc(@(p) neglog00(p,z0,z1,z2),theta_0,ops);

    lf00 = -lf00;
    lr = -2*(lf00 - lf1);
    
    disp('LR test (y2 strongly exogenous)');
    disp(['LR statistic           = ',num2str(lr) ]);
    disp(['p-value                = ',num2str(1-chi2cdf(lr,2)) ]);
    disp(' ')

    % Estimate partial model just consisting of an augmented equation 
    % for y1 by that assuming y2 is weakly exogenenous     
 
    bhat = [z1 z2 z0(:,2)]\z0(:,1);
    
    disp('Estimate of long-run parameter based on partial model')
    disp(-bhat(2)/bhat(1));

end

%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Reduced rank log-likelihood function (unrestricted)
%-------------------------------------------------------------------------
function lf = neglog1(b,z0,z1,z2)

    nobs = length(z0);
    f   = zeros(nobs,1);
    
    m1  = b(1)*(z1(:,1) - b(3)*z1(:,2)) + z2*b([4 5]);
    m2  = b(2)*(z1(:,1) - b(3)*z1(:,2)) + z2*b([6 7]);

    v1 = z0(:,1) - m1;
    v2 = z0(:,2) - m2;
    v  = [v1  v2 ];

    [k n] = size( v );
    
    omegav = v'*v/k;

    for t = 1:nobs

        f(t) = - n*0.5*log(2*pi) - 0.5*log(det(omegav)) - 0.5*v(t,:)*inv(omegav)*v(t,:)';

    end

    lf = - mean( f );
end
%-------------------------------------------------------------------------
%  Reduced rank log-likelihood function (restricted: y2 weakly exogenous)
%-------------------------------------------------------------------------
function lf = neglog0(b,z0,z1,z2)

    nobs = length(z0);
    f   = zeros(nobs,1);
    
    m1  = b(1)*(z1(:,1) - b(2)*z1(:,2)) + z2*b([3 4]);
    m2  = z2*b([5 6]);

    v1 = z0(:,1) - m1;
    v2 = z0(:,2) - m2;
    v  = [v1  v2 ];

    [k n] = size( v );
    
    omegav = v'*v/k;

    for t = 1:nobs

        f(t) = - n*0.5*log(2*pi) - 0.5*log(det(omegav)) - 0.5*v(t,:)*inv(omegav)*v(t,:)';

    end

    lf = - mean( f );
end
%-------------------------------------------------------------------------
%  Reduced rank log-likelihood function (restricted: y2 weakly exogenous)
%-------------------------------------------------------------------------
function lf = neglog00(b,z0,z1,z2)

    nobs = length(z0);
    f   = zeros(nobs,1);
    
    m1  = b(1)*(z1(:,1) - b(2)*z1(:,2)) + z2*b([3 4]);
    m2  = z2*[b(5); 0 ];

    v1 = z0(:,1) - m1;
    v2 = z0(:,2) - m2;
    v  = [v1  v2 ];

    [k n] = size( v );
    
    omegav = v'*v/k;

    for t = 1:nobs

        f(t) = - n*0.5*log(2*pi) - 0.5*log(det(omegav)) - 0.5*v(t,:)*inv(omegav)*v(t,:)';

    end

    lf = - mean( f );
end
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



