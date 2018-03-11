%=========================================================================
%
%   Program to simulate and estimate a VARMA model
%
%=========================================================================
function stsm_varma( )

    clear all
    clc

    %RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234567) );

    % Generate the data
    nobs = 1000;                 
    
    mu1 = 0.0; phi111 = 0.6; si111 = 0.2; si112= -0.5; % Equation 1    
    mu2 = 0.0; phi122 = 0.4; si121 = 0.2; si122 = 0.6; % Equation 2    

    nobs = nobs + 100;          % Allow for startup  

    y1 = zeros(nobs,1);
    y2 = zeros(nobs,1);
    v1 = randn(nobs,1);
    v2 = randn(nobs,1);


    for t = 2:nobs

        y1(t) = mu1 + phi111*y1(t-1) + v1(t) + si111*v1(t-1) + si112*v2(t-1);
        y2(t) = mu2 + phi122*y2(t-1) + v2(t) + si121*v1(t-1) + si122*v2(t-1);
    end
    y = [ trimr(y1,100,0) trimr(y2,100,0) ];     
    
    t = length(y);
    
    % Estimate the unrestricted model
    bstart = [mu1; phi111; si111; si112; mu2; phi122; si121; si122];
    options = optimset('LargeScale','off','Display','iter');                                
              
    [ theta,fvalu,~,~,~,hess ]  = ... 
                              fminunc( @(b) neglogl(b,y),bstart,options );     

    % Wald test 
    vc = (1/t)*inv( hess );
 
    r = [0   0   1   0   0   0   0   0 ;
         0   0   0   1   0   0   0   0 ;
         0   0   0   0   0   0   1   0 ;
         0   0   0   0   0   0   0   1 ] ;
    q = [ 0; 0; 0; 0 ];
    w = (r*theta - q)'*inv(r*vc*r')*(r*theta - q);

    disp( ['Wald test  = ' num2str(w) ] );
    disp( ['p-value    = ' num2str(1-chi2cdf(w,4))] );
     
    % Estimate the restricted model
    bstart = [ mu1; phi111; mu2; phi122 ];

    [ ~,fvalr,~,~,~,~ ]  = ...
                              fminunc( @(b) negloglr(b,y),bstart,options );
 
    % Likelihood Ratio test 
    fvalr = -fvalr;
    fvalu = -fvalu;
    lr = -2*(t-1)*(fvalr - fvalu);

    disp( ['LR test  = ' num2str(lr) ] );
    disp( ['p-value    = ' num2str(1-chi2cdf(lr,4))] );    
    
end
%--------------------------- Subroutines ----------------------------------
% 
%--------------------------------------------------------------------------
% Log-likelihood function for an unrestricted VARMA model
%--------------------------------------------------------------------------
function f = neglogl( b,y )

    [ t,n ] = size( y );
    e1      = zeros( t,1 );
    e2      = zeros( t,1 );
    lf      = zeros( t-1,1 );
    
    % First loop over MA part  
    for i = 2:t
        
        e1(i) = y(i,1)-b(1)-b(2)*y(i-1,1)-b(3)*e1(i-1)-b(4)*e2(i-1);
        e2(i) = y(i,2)-b(5)-b(6)*y(i-1,2)-b(7)*e1(i-1)-b(8)*e2(i-1);
    end
    e  = [ trimr( e1,1,0 ) trimr( e2,1,0 )] ;
    
    omega = e'*e/(t-1);

    for i = 1:t-1;

        lf(i) = -0.5*n*log(2*pi) - 0.5*log(det(omega)) ...
                - 0.5*e(i,:)*inv(omega)*e(i,:)';    
    end
    f = -mean( lf );
    
end
%--------------------------------------------------------------------------
% Log-likelihood function for a restricted VARMA model
%--------------------------------------------------------------------------

function f = negloglr( b,y )

    [ t,n ] = size( y );
    e1      = zeros( t,1 );
    e2      = zeros( t,1 );
    lf      = zeros( t-1,1 );
    
    % First loop over MA part  
    for i = 2:t
        
        e1(i) = y(i,1)-b(1)-b(2)*y(i-1,1);
        e2(i) = y(i,2)-b(3)-b(4)*y(i-1,2);
    end
    e  = [ trimr( e1,1,0 ) trimr( e2,1,0 )] ;
    
    omega = e'*e/(t-1);

    for i = 1:t-1;

        lf(i) = -0.5*n*log(2*pi) - 0.5*log(det(omega)) ...
                - 0.5*e(i,:)*inv(omega)*e(i,:)';    
    end
    f = -mean( lf );
    
end
