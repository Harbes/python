%=========================================================================
%
%   Program to estimate the Campbell and Shiller present value model
%
%=========================================================================

function stsm_camshiller( )

    clear all
    clc
    
    % Load data from 1933:1 to 1990:12   
    %	1. Real equity price (US)
	%	2. Real dividend (US)
 
    load campbell_shiller.mat

    % Define variables
    s = ytdata(:,1);		
    d = ytdata(:,2);
    t = length(ytdata);

    % Estimate alpha from the present value cointegrating equation
    x     = [ ones( length(d),1 ) d ];
    b_pv  = x\s;
    alpha = b_pv(2);

    % Construct variables for use in VAR
    y    = [ 100*(trimr(d,1,0) - trimr(d,0,1))  trimr(s - x*b_pv,1,0) ];
    nobs = length( y );

    % Estimate VAR(1)
    b_var = [ ones(nobs-1,1) trimr(y,0,1) ]\trimr(y,1,0);     

    % Estimate the unrestricted model
    options = optimset( 'LargeScale',           'off', ...
                        'Display',              'iter', ...
                        'MaxIter',              2000,   ...
                        'MaxFunEvals',          4000 );                                      
                 
    [ theta0,fvalu]  = fminunc( @(b) neglogl(b,y),b_var(:),options );     
      
    % Estimate the restricted model
    bstart = [ theta0(1:4);  0.95 ];

    [ theta1,fvalr]  = fminunc( @(b) negloglr(b,y,alpha),bstart,options );
 
    disp( ['Estimated discount rate       = ', num2str( 1/theta1(end)-1 ) ] );

    % Likelihood Ratio test 
    fvalr = -fvalr;
    fvalu = -fvalu;
    lr = -2*t*(fvalr - fvalu);
    dof = length(theta0) - length(theta1);

    disp( ['LR test                       = ' num2str(lr) ] );
    disp( ['Number of degrees of freedom  = ', num2str(dof) ] );
    disp( ['p-value                       = ' num2str(1 - chi2cdf(lr,4))] );

end
%--------------------------- Functions ----------------------------------
% 
%--------------------------------------------------------------------------
% Log-likelihood function for an unrestricted Campbell-Shiller model
%--------------------------------------------------------------------------
function f = neglogl( b,y )

    [ nobs,n ] = size( y );

     e  = zeros( nobs,2);    
     lf = zeros( nobs-1,1 );

     for t = 2:nobs

         e(t,1) = y(t,1) - b(1) - b(2)*y(t-1,1) - b(3)*y(t-1,2);
         e(t,2) = y(t,2) - b(4) - b(5)*y(t-1,1) - b(6)*y(t-1,2);
     end
     e = trimr(e,1,0);

     omega = e'*e/(nobs-1);

     for t = 1:nobs-1

        lf(t) = - n*0.5*log(2*pi) - 0.5*log(det(omega)) ...
                 - 0.5*e(t,:)*inv(omega)*e(t,:)';
     end
     f = -mean(lf);
end
%--------------------------------------------------------------------------
% Log-likelihood function for the restricted Campbell-Shiller model
%--------------------------------------------------------------------------
function f = negloglr( b,y,alpha )

    [ nobs,n ] = size( y );

     e  = zeros( nobs,2);    
     lf = zeros( nobs-1,1 );

     for t = 2:nobs

        e(t,1) = y(t,1) - b(1) - b(2)*y(t-1,1) - b(3)*y(t-1,2);
        e(t,2) = y(t,2) - b(4) - (-alpha*b(2))*y(t-1,1) - (1/b(5) - alpha*b(3))*y(t-1,2);

     end
     e = trimr(e,1,0);

     omega = e'*e/(nobs-1);

     for t = 1:nobs-1

        lf(t) = - n*0.5*log(2*pi) - 0.5*log(det(omega)) ...
                 - 0.5*e(t,:)*inv(omega)*e(t,:)';
     end
    f = -mean(lf);
end
