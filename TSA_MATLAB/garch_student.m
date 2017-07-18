%=========================================================================
%
%   Estimating GARCH - Student t Model
%  
%=========================================================================
function garch_student(  )

    clear all
    clc
    
    % Load data from .mat file
    
    load equity

    % Choose equity index and compute percentage return and subtract mean
    equity  = ftse;                                                 
    y       = 100*(trimr( log( equity ),1,0 ) - trimr( log( equity ),0,1 )); 
    y       = y - mean(y);
    t       = length( y );
     
    % Estimating the GARCH_t(1,1) model
    options = optimset( 'LargeScale','off',  ...
                        'Display','off',     ...
                        'MaxIter',15000,     ...
                        'MaxFunEvals',10000, ...
                        'TolFun',1e-8,       ...
                        'TolX',1e-8);
    
    start                  = [0.1 0.05 0.9 0.1 8];
    [theta,lf1,~,~,~,hess] = fminunc( @(b) lnl(b,y),start,options);

    lf1 = -lf1;
    vc  =(1/t)*inv(hess);
    
    disp('GARCH(1,1) unrestricted '); 
    disp('Estimates  Std.Errors')
    disp([theta' sqrt(diag(vc)) ]);
    disp(['Unconditional Variance = ',num2str(theta(2)/(1-theta(3)-theta(4)))]);
    disp(['Log-likelihood (unrestricted) = ',num2str(lf1) ]);
    
    % Wald test of normality: Transform the parameter so under H0 it is zero
    c  = 1/theta(5); 
    q  = 0;
    dc = -1/theta(5)^2;     % Jacobian  
    wd = t*(c - q)*inv(dc*vc(5,5)*dc)*(c - q);
    
    disp(['Wald statistic    = ',num2str(wd) ])
    disp(['p-value           = ',num2str(1-chi2cdf(wd,1)) ]);


    
     % Estimating the GARCH_t(1,1) model restricted
     start  = [0.1 0.1 8];
    [theta0,lf0,~,~,~,hess] = fminunc( @(b) lnl0(b,y),start,options);
    
    lf0 = -lf0;
    vc  = inv(hess);
    
    disp('GARCH(1,1) restricted '); 
    disp('Estimates  Std.Errors')
    disp([theta0' sqrt(diag(vc)) ]);
    disp(['Log-likelihood (unrestricted) = ',num2str(lf0) ]);
 
end

%
%--------------------------- Functions ----------------------------------
% 
%-------------------------------------------------------------------------
% Likelihood wrapper function
% This function calls the lnlt function and returns the average log-likelihood.
%-------------------------------------------------------------------------

function logl = lnl( b,y )

    logl = -mean( lnlt( b,y ) ); 
    
end

%-------------------------------------------------------------------------
% Likelihood function for a GARCH-t(1,1) model
%-------------------------------------------------------------------------

function loglt = lnlt( b,y )
    
    b = abs(b);
    t = length( y );
    u = y - b(1);  
    h = std( y )^2*ones( t,1 );
    
    for i = 2:t
        
        h(i) = b(2) + b(3)*u(i-1)^2 + b(4)*h(i-1);

    end
    loglt = - 0.5*log( h ) + log(gamma((b(5)+1)/2)/((pi*(b(5)-2)^0.5)*gamma(b(5)/2))) - ((b(5)+1)/2)*log(1+(u./sqrt( h )).^2/(b(5)-2));
    
end
%-------------------------------------------------------------------------
% Likelihood wrapper function for restricted model
%-------------------------------------------------------------------------

function logl = lnl0( b,y )

    logl = -mean( lnlt0( b,y ) ); 
    
end

%-------------------------------------------------------------------------
% Likelihood function for a GARCH-t(1,1) restricted model
%-------------------------------------------------------------------------

function loglt1 = lnlt0( b,y )
    
    b = abs(b);
    t = length( y );
    u = y - b(1);  
    h = b(2)*ones( t,1 );
    
    for i = 2:t
        
        h(i) = b(2);

    end
    loglt1 = - 0.5*log( h ) + log(gamma((b(3)+1)/2)/((pi*(b(3)-2)^0.5)*gamma(b(3)/2))) - ((b(3)+1)/2)*log(1+(u./sqrt( h )).^2/(b(3)-2));
    
end