%=========================================================================
%
%   Compute GMM estimates of the parameters of a gamma distribution
%
%=========================================================================
function gmm_gamma( )

    clear all;
    clc;
    format short;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) ); 


    t = 10;         
    % Simulate data
     y = round(5 + 2*randn(t,1));      
    % Load GAUSS data
    %y = load('table.dat','-ascii'); 

    
    % Using 2 moment conditions 
    % Zero iteration
    alpha0 = mean(y);               
    g = numgrad(@gmmcrit2,alpha0,y);           
    h = numhess(@gmmcrit2,alpha0,y);

    disp('Zero iteration');
    disp([' Parameter estimate = ', num2str(alpha0) ] );
    disp([' Objective function = ', num2str( gmmcrit2(alpha0,y) ) ] );
    disp([' Gradient           = ', num2str(g) ] );
    disp([' Hessian            = ', num2str(h) ] );
    % Newton Raphson update
    alpha1 = alpha0 - inv(h)*g;     

    
    % First iteration
    g = numgrad(@gmmcrit2,alpha1,y);           
    h = numhess(@gmmcrit2,alpha1,y);
 
    disp(' ');
    disp('First iteration');
    disp([' Parameter estimate = ', num2str(alpha1) ] );
    disp([' Objective function = ', num2str( gmmcrit2(alpha1,y) ) ] );
    disp([' Gradient           = ', num2str(g) ] );
    disp([' Hessian            = ', num2str(h) ] );
    % Newton Raphson update
    alpha2 = alpha1 - inv(h)*g;     

    % Second iteration
    g = numgrad(@gmmcrit2,alpha2,y);           
    h = numhess(@gmmcrit2,alpha2,y);
 
    disp(' ');
    disp('Second iteration');
    disp([' Parameter estimate = ', num2str(alpha2) ] );
    disp([' Objective function = ', num2str( gmmcrit2(alpha2,y) ) ] );
    disp([' Gradient           = ', num2str(g) ] );
    disp([' Hessian            = ', num2str(h) ] );
    
    v = inv(h)/t;
    disp(['Variance of alpha   = ', num2str(v) ]);
    disp(['Std error of alpha  = ', num2str(sqrt(v)) ]);

    % Iterative solution
    ops = optimset('LargeScale','off','Display','off');
    [alphahat,fc,~,~,~,H]  = fminunc(@(alpha) gmmcrit2(alpha,y),alpha0,ops);

    disp(' ');
    disp('Two moment conditions')
    disp([' Parameter estimate = ', num2str(alphahat) ] );
    disp([' Objective function = ', num2str( fc ) ] );
    v = inv(H)/t;
    disp(['Variance of alpha   = ', num2str(v) ]);
    disp(['Std error of alpha  = ', num2str(sqrt(v)) ]);
 
 
    % Using 3 moment conditions
    [alphahat,fc,~,~,~,H]  = fminunc(@(alpha) gmmcrit3(alpha,y),alpha0,ops);

    disp(' ');
    disp('Three moment conditions')
    disp([' Parameter estimate = ', num2str(alphahat) ] );
    disp([' Objective function = ', num2str( fc ) ] );
    v = inv(H)/t;
    disp(['Variance of alpha   = ', num2str(v) ]);
    disp(['Std error of alpha  = ', num2str(sqrt(v)) ]);
 
    % Estimating two-parameter gamma distribution
    theta0                 = [8 ; 2];
    [thetahat,fc,~,~,~,H]  = fminunc(@(theta) gmmcrit(theta,y),theta0,ops);

    disp(' ');
    disp('Two parameter gamma distribution')
    disp([' Estimate of alpha  = ', num2str(thetahat(1)) ] );
    disp([' Estimate of beta   = ', num2str(thetahat(2)) ] );
    disp([' Objective function = ', num2str( fc ) ] );
    v = inv(H)/t;
    disp('Covariance matrix');
    disp(v);

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
% GMM objective function - 2 moment conditions
%-------------------------------------------------------------------------
function q = gmmcrit2(alpha,y)

        m1 = y - alpha;
        m2 = y.^2 - alpha*(alpha+1);
        m  = [m1  m2];
        w  = m'*m/length(y);
        q  = 0.5*mean(m)*inv(w)*mean(m)';
end
%-------------------------------------------------------------------------
% GMM objective function - 3 moment conditions
%-------------------------------------------------------------------------
function q = gmmcrit3(alpha,y)

        m1 = y - alpha;
        m2 = y.^2 - alpha*(alpha+1);
        m3 = 1./y - 1/(alpha-1);
        m  = [m1  m2 m3];
        w  = m'*m/length(y);
        q  = 0.5*mean(m)*inv(w)*mean(m)';
end
%-------------------------------------------------------------------------
% GMM objective function - 2 parameter gamma function
%-------------------------------------------------------------------------
function q = gmmcrit(theta,y)

        alpha = theta(1);
        beta  = theta(2);
        m1    = y - alpha/beta;
        m2    = y.^2 - alpha*(alpha+1)/beta^2;
        m3    = 1./y - beta/(alpha-1);
        m     = [m1  m2  m3 ];
        w     = m'*m/length(y);
        q     = 0.5*mean(m)*inv(w)*mean(m)';
end
