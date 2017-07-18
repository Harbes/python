%=========================================================================
% 
% Program for Neyman's Smooth Goodness of Fit Test   
%   
%=========================================================================

function test_smooth( )

    clear all;
    clc;

    % Initialise random number generator    
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123457) )
      
    % Simulate data
    t = 1000;
    y = randn(t,1);        % normal (simulate under H0)
%    y = randn(t,1).^2;     % chisquared with 1 degree of freedom (simulate under H1)

    % Transform the data to uniform values and estimate my maximum likelihood
    u = normcdf(y);

    theta0 = 0.1*ones(4,1);

    % Estimate model
    options = optimset('LargeScale', 'off', 'Display', 'off');
    [theta,lnlu,~,~,~,H] = fminunc(@neglog,theta0,options,u);
    %[theta,lnlu] = fminsearch(@neglog,theta0,options,u)
    
    % LR test
    lnl1 = -lnlu;         
    lnl0 = 0.0;     % Under the null lnl is zero as f under H0 is 1 so lnl=0
    lr   = -2*t*(lnl0 - lnl1);
    dof = 4;
    disp(['LR test statistic    = ' num2str(lr) ]);
    disp(['Degrees of freedom   = ' num2str(dof) ]);
    disp(['P-value              = ' num2str(1-chi2cdf(lr,dof)) ]);
    disp(' ');


    % Wald test
	wd = t*theta'*(inv(H))*theta;
    disp(['Wald statistic            = ' num2str(wd) ]);
    disp(['P-value                   = ' num2str(1-chi2cdf(wd,dof)) ]);
    disp(' ')

    % LM test
    z = u - 0.5;
    phi1 = sqrt(3)*2*z;
    phi2 = sqrt(5)*(6*z.^2 - 0.5);
    phi3 = sqrt(7)*(20*z.^3 - 3*z);
    phi4 = 3*(70*z.^4 - 15*z.^2 + 3/8);

    lm = (sum(phi1)/sqrt(t))^2 + (sum(phi2)/sqrt(t))^2 ...
        + (sum(phi3)/sqrt(t))^2 + (sum(phi4)/sqrt(t))^2;
 
    disp(['LM statistic                = ' num2str(lm) ]);
    disp(['P-value                     = ' num2str(1-chi2cdf(lm,dof)) ]);
    
end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
% Program used to compute normalizing constant numerically
%-------------------------------------------------------------------------
function y = normconst(u,theta)

    z = u - 0.5;
    t = theta(1)*sqrt(3)*2*z + theta(2)*sqrt(5)*(6*z.^2 - 0.5) ...
        + theta(3)*sqrt(7)*(20*z.^3 - 3*z) + theta(4)*3*(70*z.^4 - 15*z.^2 + 3/8);
    y = exp(1 + t);
end
%-------------------------------------------------------------------------
% Log-likelihood function
%-------------------------------------------------------------------------
function lf =  neglog(b,u)

    z = u - 0.5;
    c  = quad(@(u) normconst(u,b),0,1);            
    f  = 1 + b(1)*sqrt(3)*2*z + b(2)*sqrt(5)*(6*z.^2 - 0.5) ...
        + b(3)*sqrt(7)*(20*z.^3 - 3*z) + b(4)*3*(70*z.^4 - 15*z.^2 + 3/8) - log(c);
    lf = -mean( f );
    
end


