%=========================================================================
%
%   Program to estimate a duration model of strikes
%   using the Kennan (1985) data set.
%
%=========================================================================
function discrete_strike( )
    
    clear all
    clc

    % Load data: US strike data from 1968 to 1976
    load strike.mat

    y   = data(:,1);    % Duration of strike in days    
    out = data(:,2);    % Unanticipated output     
    t   = length(y);
    x   = [ ones(t,1)   out ];

    % Estimate the Weibull model by MLE 
    ops = optimset('LargeScale','off','Display','off');

    theta_0       = [ x\log(y) ; 1];
    [ theta1,l1 ] = fminsearch(@(b) neglog(b,y,x),theta_0,ops);
    
    l1 = -l1;   
    disp(['Log-likelihood (unrestricted) = ',num2str(l1) ]);
    disp(['(T-1)xLog-likelihood        = ',num2str((t-1)*l1) ]);
    disp(' ');

    % Estimate the Exponential model by MLE
    theta_0       = [ x\log(y) ];
    [ theta0,l0 ] = fminsearch(@(b) neglog0(b,y,x),theta_0,ops);

    l0 = -l0;   
    disp(['Log-likelihood (restricted) = ',num2str(l0) ]);
    disp(['(T-1)xLog-likelihood        = ',num2str((t-1)*l0) ]);
    disp(' ');

    % Likelihood ratio test 

    lr = -2*t*(l0 - l1);

	disp(['LR Statistic         = ',num2str(lr)]);
    disp(['p-value              = ',num2str(1-chi2cdf(lr,size(x,2)-1))]);
    disp(' ');

end
%
%--------------------------- Functions -----------------------------------
% 

%-------------------------------------------------------------------------
%  Unrestricted log-likelihood function of duration model
%-------------------------------------------------------------------------
function lf = neglog(b,y,x)

	m  = x*b(1:2);   
    s2 = abs(b(3));   
    z  = (log(y) - m)./sqrt(s2);
    lf = -mean( 0.5*log(s2) + z - exp(z) );

end

%-------------------------------------------------------------------------
%  Restricted log-likelihood function of duration model
%-------------------------------------------------------------------------
function lf = neglog0(b,y,x)

	m  = x*b(1:2);   
    s2 = 1;   
    z  = (log(y) - m)./sqrt(s2);
    lf = -mean( 0.5*log(s2) + z - exp(z) );

end

