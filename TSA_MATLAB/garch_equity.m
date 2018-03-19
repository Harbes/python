%=========================================================================
%
%   Estimating GARCH Models by ML
%
%=========================================================================

function garch_estimate(  )
%%
    clear all
    clc
   %% 
    % Load data
    load equity
%%
    % Choose equity index and compute percentage return and subtract mean
    equity  = ftse; % ftse, dow, nikkei
    
    y       = 100*(trimr( log( equity ),1,0 ) - trimr( log( equity ),0,1 )); 
    y       = y - mean(y);
    t       = length( y );
     
    % Estimating the GARCH(1,1) model
    ops = optimset( 'LargeScale','off','Display','off');
    start  = [ 0.05 0.1 0.9 ];
    %%
    tic
    [theta,lf1,~,~,~,hess] = fminunc( @(b) neglog(b,y),start,ops);
    toc
%%
    clear start
    lf1 = -lf1;
    
    disp(' ');
    disp('GARCH(1,1) Results');
    disp(['alpha_0                                 = ',num2str(theta(1)) ]);
    disp(['alpha_1                                 = ',num2str(theta(2)) ]);
    disp(['beta_1                                  = ',num2str(theta(3)) ]);
    disp(['Log-likelihood function (unrestricted)  = ',num2str(lf1) ]);    
    
    % Compute conditional variance at optimal parameters
    h = std( y )^2*ones( t,1 );
    u = y - theta(1);
    
    for i = 2:t
        
        h(i) = theta(1) + theta(2)*u(i-1)^2 + theta(3)*h(i-1);

    end

    %*********************************************************************
    %**     Generate graph of conditional variance
    %*********************************************************************
    figure(1);
    clf;

    plot( 1:length(y), h,'-k'); 
    ylabel('Sigma^2');
    xlabel('t');

    % Comparison of standard errors 
    ih  = inv(hess);
    seh = sqrt(diag((1/t)*ih));  
    g   = numgrad( @lnlt,theta',y);
    j   = g'*g/t;
    sej = sqrt(diag( (1/t)*inv(j)));
    seq = sqrt(diag( (1/t)*( ih*j*ih)));

    disp(' ');
    disp( 'Comparison of Standard Errors' );
    disp( '    Hessian   OPG       QMLE  ');
    disp( [ seh sej seq ] );

    % Estimating the restricted ARCH(1) model
    start        = [ 0.05 0.1  ];
    [theta0,lf0] = fminunc( @(b) neglog0(b,y),start,ops);

    clear start
    lf0 = -lf0;
    disp(' ');
    disp('ARCH(1) Results');
    disp(['alpha_0                                 = ',num2str(theta0(1)) ]);
    disp(['alpha_1                                 = ',num2str(theta0(2)) ]);
    disp(['Log-likelihood function (unrestricted)  = ',num2str(lf0) ]);    

    % Likelihood Ratio test    
    lr  = -2*t*(lf0 - lf1);
    pv = 1 - chi2cdf(lr,1);   
    disp(['LR  test   = ',num2str(lr) ]);
    disp(['p-value    = ',num2str(pv) ])

    % Estimating the IGARCH;
    start         = [ 0.05 0.1]; 
    [thetai,lfi]  = fminunc( @(b) neglogi(b,y),start,ops);
    
    thetai = abs(thetai);
    lfi = -lfi;
    disp(' ');
    disp('GARCH(1,1) Results');
    disp(['alpha_0                           = ',num2str(thetai(1)) ]);
    disp(['alpha_1                           = ',num2str(thetai(2)) ]);
    disp(['beta_0                            = ',num2str(1-thetai(2)) ]);
    disp(['Log-likelihood function (IGARCH)  = ',num2str(lfi) ]);    

end
%
%--------------------------- Functions  ----------------------------------
% 
%-------------------------------------------------------------------------
% Wrapper function for an GARCH model
%-------------------------------------------------------------------------

function logl = neglog( b,y )

    logl = -mean( lnlt( b,y ) ); 
    
end

%-------------------------------------------------------------------------
% Likelihood function for a GARCH(1,1) model
%-------------------------------------------------------------------------

function loglt = lnlt( b,y )
    
    b = abs(b);
    t = length( y );
    u = y ;  
    h = std( y )^2*ones( t,1 );
    
    for i = 2:t
        
        h(i) = b(1) + b(2)*u(i-1)^2 + b(3)*h(i-1);

    end
    loglt = - 0.5*log( 2*pi ) - 0.5*log( h ) - 0.5*(u./sqrt( h )).^2;
    
end

%-------------------------------------------------------------------------
% Wrapper function for an ARCH model
%-------------------------------------------------------------------------
function logl = neglog0( b,y )

    logl = -mean( lnlt0( b,y ) ); 
    
end
%-------------------------------------------------------------------------
% Likelihood function for a ARCH(1) model
%-------------------------------------------------------------------------
function loglt = lnlt0( b,y ) 
 
    b = abs(b);
    t = length( y );
    u = y;  
    h = std( y )^2*ones( t,1 );
    
    for i = 2:t
        
        h(i) = b(1) + b(2)*u(i-1)^2;

    end
    loglt = - 0.5*log( 2*pi ) - 0.5*log( h ) - 0.5*(u./sqrt( h )).^2;
 
end

%------------------------------------------------------------------------
% Wrapper function for an IGARCH model
%-------------------------------------------------------------------------
function logl = neglogi( b,y )

    logl = -mean( lnlti( b,y ) ); 
    
end
%-------------------------------------------------------------------------
% Likelihood function for a IGARCH* model
%-------------------------------------------------------------------------
function loglt = lnlti( b,y )
    
    b = abs(b);
    t = length( y );
    u = y ;  
    h = std( y )^2*ones( t,1 );
    
    for i = 2:t
        
        h(i) = b(1) + b(2)*u(i-1)^2 + (1-b(2))*h(i-1);

    end
    loglt = - 0.5*log( 2*pi ) - 0.5*log( h ) - 0.5*(u./sqrt( h )).^2;
    
end

