%=========================================================================
%
%   Program to perform tests of the stationary distribution of the 
%   interest rate based on the gamma distribution
%
%=========================================================================
function test_interest(  )

    clear all
    clc
    
    % Load data (5505x4 array called eurodata, 1 Jun 1973 - 25 Feb 1995)
    %   1. year
    %   2. day
    %   3. date stamp
    %   4. interest rates
    load eurodollar.mat

    r = eurodata(:,4);
    t = length(r);

    % Estimate the model
    start = [1 ; 1];
    ops   = optimset('LargeScale', 'off', 'Display', 'off'); 
    
    [bhat,~,~,~,~,hess] = fminunc(@(b) neglog(b,r),start,ops);
     
    vc    = (1/t)*inv(hess);
	nu    = bhat(1);
    omega = bhat(2);
    
    disp( 'Parameter estimates')
    disp( ['nu           = ',num2str(nu) ]);
    disp( ['omega        = ',num2str(omega) ]);
    
    disp(' ')
    disp( 'Hessian matrix')
    disp( hess );

    % --------------------------------------------------------------
    % Estimate the mean and its standard error    
    mu = nu/omega;
    d  = [(1/omega)   (-nu/omega^2)];
    vr = d*vc*d';
    se = sqrt(vr);

    disp(' ')
    disp(['Mean              = ',num2str(mu) ]);     
    disp(['Std Error of Mean = ',num2str(se) ]);
    
    % Wald test of the mean 
    c = mu;
    q = 0.1;
    wd = t*(c - q)'*inv(d*(inv(hess))*d')*(c - q);

    disp(' ')
	disp(['Wald statistic = ',num2str(wd) ]);
    disp(['p-value        = ',num2str(1-chi2cdf(wd,1)) ]);
    
    % --------------------------------------------------------------
    % Estimate the variance and its standard error   
    s2 = nu/omega^2;
    d  = [(1/omega^2)   (-2*nu/omega^3)];
    vr = d*vc*d';
    se = sqrt(vr);

    disp(' ')
    disp(['Variance              = ',num2str(s2) ]);     
    disp(['Std Error of Variance = ',num2str(se) ]);

    
    % Wald test of the variance   
    c = s2;
    q = 0.001;
    wd = t*(c - q)'*inv(d*(inv(hess))*d')*(c - q);

    disp(' ')
    disp(['Wald statistic = ',num2str(wd) ]);
    disp(['p-value        = ',num2str(1-chi2cdf(wd,1)) ]);

    % --------------------------------------------------------------
    % Estimate skewness and its standard error    
    k3   = 2/sqrt(nu);
    d    = [(-1/nu^(3/2))  0];
    vr = d*vc*d';
    se = sqrt(vr);

    disp(' ')
    disp(['Skewness              = ',num2str(k3) ]);     
    disp(['Std Error of skewness = ',num2str(se) ]);

    % Wald test of the skewness   
    c = k3;
    q = 0.0;
    wd = t*(c - q)'*inv(d*(inv(hess))*d')*(c - q);
    
    disp(' ')
    disp(['Wald statistic = ',num2str(wd) ]);
    disp(['p-value        = ',num2str(1-chi2cdf(wd,1)) ]);

    % --------------------------------------------------------------
    % Estimate kurtosis and its standard error    
    k4   = 3 + 6/nu;
    d    = [-6/nu^2  0];
    vr = d*vc*d';
    se = sqrt(vr);

    disp(' ')
    disp(['Kurtosis              = ',num2str(k4) ]);     
    disp(['Std Error of kurtosis = ',num2str(se) ]);

    % Wald test of the kurtosis   
    c = k4;
    q = 3;
    wd = t*(c - q)'*inv(d*(inv(hess))*d')*(c - q);
    
    disp(' ')
    disp(['Wald statistic = ',num2str(wd) ]);
    disp(['p-value        = ',num2str(1-chi2cdf(wd,1)) ]);
   
end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
% Likelihood function for stationary distribution of CIR model
%-------------------------------------------------------------------------
function lf = neglog(b,y)

    nu = abs(b(1));
    om = abs(b(2));
    f  = nu*log( om ) - gammaln( nu ) + (nu-1)*log(y) - om*y;  
    lf = -mean(f);
end


