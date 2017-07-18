%=========================================================================
%
%   Program to estimate a threshold autoregressive model based on 
%   the stationary and transitional distributions
%
%=========================================================================
function max_tar( )

    clear all
    clc

	RandStream.setDefaultStream( RandStream('mt19937ar','seed',12357) )

    t     = 250; 
    theta = 0.5;  

    % Simulate data 
    y = zeros(t,1);

    for i=2:t

        y(i) = theta*abs(y(i-1)) + randn;

    end

    % Estimate by least squares       
    xvar = trimr(abs(y),0,1);
    yvar = trimr(y,1,0); 
    bols = xvar\yvar;

    % Estimate by Ml applied to the stationary distribution       
    start = bols;
    ops   = optimset('LargeScale','off','Display','off');
    [ bhat,~,~,~,~,hess] = fminunc(@(b) lnlstat( b,y ),start,ops );


    disp('True value of theta')
    disp( theta )
    
    disp('OLS Results')
    disp( bols );
    
    disp('Stationary distribution results')
    disp( bhat )


    % Estimate by Ml applied to the transtional distribution       
    start = bols;
    ops   = optimset('LargeScale','off','Display','off');
    [ bhat,~,~,~,~,hess] = fminunc(@(b) lnltrans( b,y ),start,ops );

    disp('Transitional distribution results')
    disp( bhat )
    
end
% ------------------------ Functions ------------------------------------%
%
%-------------------------------------------------------------------------
%  Log-likelihood of a tar model: transitional distribution
%-------------------------------------------------------------------------
function lf = lnltrans(b,y)

    m = b*abs(trimr(y,0,1));
    s = 1;
    z = ( trimr(y,1,0) - m ) ./ s;
	f = log( normpdf(z)/s );                                  

    lf = -mean(f);    
end
%-------------------------------------------------------------------------
%  Log-likelihood of a tar model: stationary distribution
%-------------------------------------------------------------------------
function lf = lnlstat(b,y)

        m = 0.0;
        s = 1/sqrt(1 - b^2);
        z = ( y - m ) ./ s;
        f = 2*normpdf(z)*(1/s).*normcdf(b*y);    
        
        lf = -mean( log(f) );
end

