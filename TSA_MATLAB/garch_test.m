%=========================================================================
%
%      Program to test for and estimate an ARCH(p) model
%
%=========================================================================
function garch_test( )

    clear all;
    clc;

    % Load data
    load equity;
 
    % Choose index
    y = ftse;   % ftse, dow, nikkei
    y = 100*(log(trimr(y,1,0)) - log(trimr(y,0,1)));
    y=y-mean(y);

    % Set ARCH order
    p = 2;
    
    % Test for ARCH(p)
    [lm,pv]=testarch(y,p);
    
    % Estimate ARCH(p)
    ops        = optimset('LargeScale','off','Display','off');
    start      = rand(p+1,1);
    [theta,lf] = fminunc( @(b) neglog(b,y,p),start,ops);
    theta      = abs(theta);    
    
    disp('ARCH Estimation and Testing');
    disp('---------------------------')
    disp(['ARCH order          = ',num2str(p) ]);
    disp(['alpha_0             = ',num2str(theta(1)) ]);
    disp(['alpha_1             = ',num2str(theta(2)) ]);
    disp(['Likelihood function = ',num2str(-lf) ]);
    disp(' ');
    disp(['LM test             = ',num2str(lm) ]);
    disp(['p-value             = ',num2str(pv) ]);
    
end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Test for ARCH(maxlag)
%-------------------------------------------------------------------------
function [ lm,pv ] = testarch(data,maxlag)

    % Create lags of y
    tmp = lagmatrix(data,1:maxlag);
    tmp(any(isnan(tmp),2),:) = [];

    y = trimr(data,maxlag,0).^2;
    x = [ ones(length(tmp),1) tmp.^2 ];
    b = x\y;
    e = y - x*b;

    t   = length(y);    
    yc  = y-mean(y);
    rss = e'*e;
    tss = yc'*yc;
    r2  = 1- rss/tss;
    lm  = t*r2;
    pv  = 1-chi2cdf(lm,maxlag);
end
%-------------------------------------------------------------------------
% Likelihood function for an ARCH(p) model
%-------------------------------------------------------------------------
function lf = neglog( b,data,maxlag )

    b = abs(b);
    
    % Compute u and lags of u
    u   = data; 
    tmp = lagmatrix(data,1:maxlag);
    tmp(any(isnan(tmp),2),:) = [];
    
    % Compute h
    h = [ ones(length(tmp),1) tmp.^2 ]*b;
    h = [std(u)^2*ones(maxlag,1) ; h];

    % Likelihood function
    f  = - 0.5*log( 2*pi ) - 0.5*log( h ) - 0.5*(u./sqrt( h )).^2;
    lf = -mean( f );
    
end
