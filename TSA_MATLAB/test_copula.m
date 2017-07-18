%=========================================================================
%
%   Program to estimate a bivariate gaussian copula for asset returns
%   Asset price data from 6 August 2010 to 2 January 2001 (note that the
%   data are in reverse order ie from recent to past)
%
%=========================================================================
function test_copula( )

    clear all
    clc

    % Load data
    load diversify.mat

    t = 2413;
    
    % Select appropriate sample  
    pt_apple     = pt_apple(1:t);
    pt_ford      = pt_ford(1:t);

    % Compute percentage returns  
    r_apple = 100*diff(log(pt_apple)); 
    r_ford  = 100*diff(log(pt_ford));

    y = [r_apple r_ford];
    t =length(y);

    % Compute statistics
    m = mean(y);
    s = std(y);
    c = corrcoef(y);
    r = c(1,2);

    % Estimate parameters of the copula 
    start = [ m' ; s' ; r ];
    ops   = optimset('LargeScale','off','Display','off');
    
    [ bhat,~,~,~,~,hess] = fminunc(@(b) neglog( b,y ),start,ops );

    vc = (1/t)*inv(hess);
    
    % Wald test of independence        
    wd = (bhat(5) - 0)^2/vc(5,5);
    disp(' ');
    disp( ['Wald statistic     = ',num2str(wd) ]);
    disp( ['P-value            = ',num2str(1-chi2cdf(wd,1)) ]);

end

%
% ------------------------ Functions ------------------------------------%
%
%-------------------------------------------------------------------------
%  Copula log-likelihood function
%-------------------------------------------------------------------------
function lf = neglog(b,y)

        m1 = b(1);                 % Asset 1
        s1 = abs(b(3));
        z1 = (y(:,1) - m1)/s1;
        f1 = normpdf(z1)/s1;

        m2 = b(2);                 % Asset 2 
        s2 = abs(b(4));
        z2 = (y(:,2) - m2)/s2;
        f2 = normpdf(z2)/s2;

        r  = b(5);                 % Dependence

        lt = -log(2*pi) - log(s1*s2) - 0.5*log(1 - r^2) ... 
            - 0.5*(z1.^2 - 2*r*z1.*z2 + z2.^2)/(1 - r^2) + log(f1) + log(f2);

        lf = -mean(lt);
        
end