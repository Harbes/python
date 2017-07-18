%==========================================================================
%
%   	Program to estimate a DECO model of US yields
%
%==========================================================================
function mgarch_deco( )

    clear all
    clc

    global RBAR
    
    
    % Load data: US daily yields 3-Jan-2000 to 21-Aug-2006
    load daily_finance
   
    % Data manipulation
    data = yields(:,[1 5 10 15 20]);
    y    = 100*(trimr(data,1,0) - trimr(data,0,1));
    y    = bsxfun(@minus, y, mean(y));
    
    [ t,n ] = size( y );  
    
    RBAR = zeros( t,1 );
    
    % Estimate univariate GARCH models for each variable   
    ops   = optimset( 'LargeScale','off','Display','off' );
    b1 = zeros( 3,n );
    h1 = zeros( t,n );
    h  = zeros( t,1 );

    for k = 1:n 
        
        bstart      = [0.1  -2  2 ];
        b1(:,k)     = fminunc( @(b) negloglgarch(b,y(:,k)),bstart,ops);   
        %b1(:,k)     =[ tmp(1) normcdf(tmp(2)) normcdf(tmp(3)) ];
        [~,h1(:,k)] = negloglgarch( b1(:,k),y(:,k) );

    end
    b1([2 3],:) = normcdf( b1([2 3],:) );
    disp( 'Garch parameters' );
    disp( b1 );
    
    % Estimate the DECO model
    bstart = [ -2 2  ];
    bc = fminunc( @(b) neglogldeco(b,y,h1),bstart,ops);    
    bc = normcdf( bc );
    disp( bc );
    
    start =           [     0.999314709903659
                            0.152115338735339
                            0.843167761478777
                            0.789746885859107
                            0.064814542202766
                            0.912043327137570
                            0.630618257788049
                            0.054182853876372
                            0.920133291548674
                            0.815035172227143
                            0.055163278439432
                            0.912538859743870
                            0.796164122428437
                            0.055205181793798
                            0.911644837026472
                            0.180319393683112
                            0.953053140443055];
    
    ops =     optimset( 'LargeScale','off',  ...
                        'Display','iter',     ...
                        'MaxIter',15000,     ...
                        'MaxFunEvals',10000, ...
                        'TolFun',1e-5,       ...
                        'TolX',1e-7);

    bf = fminunc( @(b) neglogfull(b,y),start,ops);
    
    ind = [ 2 3 5 6 8 9 11 12 14 15 16 17 ];
    bf(ind) = normcdf(bf(ind));
    disp([[b1(:); bc(:)] bf]);
    
    
%     %*********************************************************************
%     %**     Generate graph of conditional variances
%       This was set up to read from the excel file so is commented out
%     %*********************************************************************
%     % Get daily dates from excel file
%     dvec     = datenum( txt(:,1), 'mm/dd/yyyy' );
%     vec      = dvec(2:end);
%     
%     % Switch off TeX interpreter and clear figure
%     set(0,'defaulttextinterpreter','none');
%     figure(1);
%     clf;
%     
%     plot( vec,RBAR,'-k','LineWidth',0.75);
%     datetick('x','yyyy'); 
%     %title(' (a) Variance of 3-month rate' );
%     ylabel('Average Correlation');
%     xlabel('$t$');
%     set(gca,'YLimMode','manual')
%     set(gca,'Ylim',[0 1]);
%     set(gca,'YTick',[ 0.2 0.4 0.6 0.8 1] );
%     box off
% 
%     
%     
%     %laprint(1,'deco','options','factory');
%     
    
end

%
%--------------------------- Functions ----------------------------------
% 
%-------------------------------------------------------------------------
% Likelihood wrapper function
% This function calls the lnlt function and returns the average log-likelihood.
%-------------------------------------------------------------------------

function [lf,h] = negloglgarch( b,y ) 
    
    u = y ;    
    h = recserar(b(1) + normcdf(b(2))*trimr([0.0;u.^2],0,1),std(u)^2,normcdf(b(3)));  
    f = - 0.5*log( 2*pi ) - 0.5*log( h ) - 0.5*(u./sqrt( h )).^2;
    
    lf = -mean( f );
    
end

%-------------------------------------------------------------------------
% Likelihood wrapper function
% This function calls the lnlt function and returns the average log-likelihood.
%-------------------------------------------------------------------------
function lf = neglogldeco( b,y,hv )

    b       = normcdf( b );    
    [ t,n ] = size( y );
    f       = zeros( t,1 );
    
    u    = y;              % Residuals                                                        
    z    = u./sqrt(hv);    % Standardised residuals                                          
    qbar = cov( z );
    q    = qbar;                     

    for i = 1:t
        
        % Diagonal matrix of conditional standard deviations
        s  = diag( sqrt(hv(i,:)) ); 
        
        % Conditional correlation matrix 
        tmp = inv( diag( sqrt(diag(q)) ) );
        r   = tmp*q*tmp;     
        
        % Compute mean of covariances
        ind  = vec( tril( r,-1 ) ) ~= 0 ;
        rbar = mean( r(ind) );                     
        
        % Redefine r for DECO    
        r = (1 - rbar)*eye(n) + rbar*ones(n,n);  
        
        % Determinant and inverse of r for DECO
        rdet = (1 - rbar)^(n-1) * ( 1 + (n - 1)*rbar);                                              
        rinv = (1/(1 - rbar))*eye(n) - (rbar/( (1-rbar)*(1 + (n - 1)*rbar) ) )*ones(n,n);       
               
        % Log of the likelihood at observation i   
        f(i) = - 0.5*log(rdet) - 0.5*z(i,:)*rinv*z(i,:)';% + 0.5*z(i,:)*z(i,:)';   
        
        % Update Q
        q  = abs(1 - b(1) - b(2))*qbar + b(1)*z(i,:)'*z(i,:) + b(2)*q;   
                 
    end
    lf = -mean( f );
end

%-------------------------------------------------------------------------
% Likelihood wrapper function for FULL DECO model
% This function calls the lnlt function and returns the average log-likelihood.
%-------------------------------------------------------------------------
function lf = neglogfull( b,y )

    ind = [ 2 3 5 6 8 9 11 12 14 15 16 17 ];
    b(ind) = normcdf(b(ind));
    
    global RBAR

    [ t,n ] = size( y );                                                            
    f      = zeros( t,1 );                                                       
    u = y;                                                                               
    % Construct conditional variances
    hv(:,1) = recserar(b(1)+b(2)*trimr([0.0;u(:,1).^2],0,1),std(u(:,1))^2,b(3));
    hv(:,2) = recserar(b(4)+b(5)*trimr([0.0;u(:,2).^2],0,1),std(u(:,1))^2,b(6));      
    hv(:,3) = recserar(b(7)+b(8)*trimr([0.0;u(:,3).^2],0,1),std(u(:,1))^2,b(9));
    hv(:,4) = recserar(b(10)+b(11)*trimr([0.0;u(:,4).^2],0,1),std(u(:,1))^2,b(12));     
    hv(:,5) = recserar(b(13)+b(14)*trimr([0.0;u(:,5).^2],0,1),std(u(:,1))^2,b(15));    
    
    % Unconditional covariance matrix  standardized residuals
    z    = u./sqrt( hv );                                                
    qbar = cov( z );
    q    = qbar;                                                                        
    
    for i = 1:t
        
        % Diagonal matrix of conditional standard deviations
        s  = diag( sqrt(hv(i,:)) ); 
        
        % Conditional correlation matrix 
        tmp = inv( diag( sqrt(diag(q)) ) );
        r   = tmp*q*tmp;     
        
        % Compute mean of covariances
        ind  = vec( tril( r,-1 ) ) ~= 0 ;
        rbar = mean( r(ind) );     
        RBAR(i) = rbar;
        
        % Redefine r for DECO    
        r = (1 - rbar)*eye(n) + rbar*ones(n,n);  
                   
        % Update conditional variance-covariance matrix  
        h  = s*r*s;                                                                                    
       
        % Likelihood function
        f(i) = - 0.5*n*log(2*pi) - 0.5*log(det(h)) - 0.5*u(i,:)*inv(h)*u(i,:)'; 
       
        % Update Q
        q  = abs(1 - b(16) - b(17))*qbar + b(16)*z(i,:)'*z(i,:) + b(17)*q;   
    end
    lf = -mean( f );
      
end

%--------------------------------------------------------------------------
% Column stack a matrix to agree with GAUSS VEC command
%--------------------------------------------------------------------------
function y = vec(x)

    y = x(:);

end




