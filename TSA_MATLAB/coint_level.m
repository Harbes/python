%=========================================================================
%
%   Program to estimate a bivariate vecm of the 
%    term structure of interest rates with level effects on the variance
%=========================================================================
function coint_level( )

    clear all
    clc
    
    % Load data and choose 10-year and 1-year yields
    load usmacro.mat   
    yield = [r10yr  r1yr];

    z0   = trimr(yield,1,0) - trimr(yield,0,1);
    z1   = trimr(yield,0,1);

    t    = length(z0);
    nobs = length(yield);

    % Starting values
    b = [ones(nobs,1) yield(:,2)]\yield(:,1);
    u = yield(:,1) - [ ones(nobs,1) yield(:,2)]*b;
    a = trimr(u,0,1)\z0;   
    v = z0 - trimr(u,0,1)*a;

    % Estimate levels effect model
    start = [ b ; a' ; vech(chol(cov(v))') ; 0.0 ; 0.0 ];
    ops    = optimset('LargeScale', 'off', 'Display', 'iter');     
    [theta,logl,~,~,~,hess] = fminunc(@(p) neglog(p,z0,z1),start,ops);


    disp('Levels effect model estimates')        
    disp('-----------------------------')
    disp(['Log likelihood function = ',num2str(-logl) ])
    disp(['beta_c         = ',num2str(theta(1)) ])
    disp(['beta_r         = ',num2str(theta(2)) ])
    disp(['alpha_1        = ',num2str(theta(3)) ])
    disp(['alpha_2        = ',num2str(theta(4)) ])

    
    vc = (1/t)*inv(hess);
 
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


    % Wald test of levels effect
    r  = [ 0  0  0  0  0  0  0  1  0 ;
           0  0  0  0  0  0  0  0  1 ];
    q  = [ 0 ; 0] ;
    wd = (r*theta - q)'*inv(r*vc*r')*(r*theta - q);

    disp('Wald test of levels effect');
    disp(['Wald statistic           = ',num2str(wd) ]);
    disp(['p-value                  = ',num2str(1-chi2cdf(wd,2)) ]);
    disp(' ')

end

%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Reduced rank log-likelihood function: unrestricted
%-------------------------------------------------------------------------
function lf = neglog(b,z0,z1)

    nobs = length(z0);
    f   = zeros(nobs,1);

    v1 = z0(:,1) - b(3)*(z1(:,1) - b(1) - b(2)*z1(:,2));
    v2 = z0(:,2) - b(4)*(z1(:,1) - b(1) - b(2)*z1(:,2));
    v  = [v1  v2];
    

    sd  =  [ b(5)    0      ;          
             b(6)    b(7) ] ;       

    % Level effects parameters     
    kappa  =  [ b(8)    b(9)] ;       

    [~,n] = size( v );

    for t = 1:nobs
        
        %  Specify level-effects covariance structure using Choleski          
        l      = bsxfun(@times,sd,(z1(t,:).^kappa));   
        omegav = l*l'; 


        f(t) = - n*0.5*log(2*pi) - 0.5*log(det(omegav)) - 0.5*v(t,:)*inv(omegav)*v(t,:)';

    end

    lf = - mean( f );
end
