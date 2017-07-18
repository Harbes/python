%=========================================================================
%
%  Program to compare the efficiency properties of MLE and OLS
%  in the AR(1) regression model in the case where the
%  explanatory variable is a constant. In this case the OLS estimator is
%  the sample mean
%
%=========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) )
rho = 0.6;
sig2 = 10;
T    = [5, 50, 500];

%   Compute OLS variance
for k = 1:length(T) 
    t = T(k);
    sum = 0.0;
    for i = 1:t-1
        for j = 1:t-i
            sum = sum + rho^i;
        end
    end
    var_ols = sig2*(t + 2*sum)/((1-rho^2)*t^2);

    %  OLS estimator asymptotic variance based on Harvey (1990, p.197)
    sum_h = (t*(1-rho^2) - 2*rho*(1-rho^t))/((1-rho)^2);
    var_ols_harvey = sig2*(sum_h)/( (1-rho^2)*t^2 );

    % Compute MLE variance
    var_mle = sig2/((t-1)*(1-rho)^2);

    %  Print results
    disp(['Sample size              = ', num2str(T(k)) ]);
    disp(['Variance of OLS          = ', num2str(var_ols) ]);
    disp(['Variance of OLS (Harvey) = ', num2str(var_ols) ]);
    disp(['Variance of MLE          = ', num2str(var_mle) ]);
    disp(['Efficiency (OLS/MLE)     = ', num2str(var_ols/var_mle) ]);

end

