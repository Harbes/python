%=========================================================================
%
%      Program to estimate level effect in interest rates by GMM
%
%=========================================================================

function gmm_level( )
    
    clear all;
    clc;


    % Load data --- monthly December 1946 to February 1991                                          
    %     0,1,3,6,9 mth and 10 yrs                                        
    load level.mat;
    xdata = level;

    %     Choose the interest rate series         
    rt  = xdata(:,3);
    drt = trimr(rt,1,0) - trimr(rt,0,1);
    r1t = trimr(rt,0,1);
    t   = length(drt);

    lmax = 5;    %  Length of the Newey-West Covariance matrix estimation                            

    % Estimate the model to obtain consistent estimates
    flag = 0;
    ops  = optimset('LargeScale','off','Display','off');
    b0   = [0.1;0.1;0.1;1.0];              
    bgmm = fminunc(@(b) q(b,drt,r1t,lmax,flag),b0,ops);

    % Now estimate with optimal weighting matrix
    flag = 1;
    bgmm = fminunc(@(b) q(b,drt,r1t,lmax,flag),bgmm,ops);

    % Compute optimal weigthing matrix at GMM estimates   
    [obj,w] = q(bgmm,drt,r1t,lmax,flag);

    % Compute standard errors of GMM estimates
    dg = numgrad(@meaneqn,bgmm,drt,r1t);
    v  = dg'*inv(w)*dg;
    cov = inv(v)/t;
    se = sqrt(diag(cov));

    disp(' ');
	disp(['The value of the objective function  = ', num2str(obj) ]);
    disp(['J-test                               = ', num2str(t*obj) ]);
    disp('Estimates     Std err.   t-stats');
    disp( [ bgmm  se  bgmm./se ])
    disp(['Newey-West estimator with max lag    = ', num2str(lmax) ]);
    disp(' ');

    % Test of gam = 0.0,0.5,1.0,1.5 
    stat = (bgmm(4) - 0.0)^2/cov(4,4);
    disp(['Test of (gam=0.0) = ', num2str(stat) ]);
    disp(['p-value           = ', num2str(1-chi2cdf(stat,1)) ]);
    disp(' ');


    stat = (bgmm(4) - 0.5)^2/cov(4,4);
    disp(['Test of (gam=0.5) = ', num2str(stat) ]);
    disp(['p-value           = ', num2str(1-chi2cdf(stat,1)) ]);
    disp(' ');

    stat = (bgmm(4) - 1.0)^2/cov(4,4);
    disp(['Test of (gam=1.0) = ', num2str(stat) ]);
    disp(['p-value           = ', num2str(1-chi2cdf(stat,1)) ]);
    disp(' ');

    stat = (bgmm(4) - 1.5)^2/cov(4,4);
    disp(['Test of (gam=1.5) = ', num2str(stat) ]);
    disp(['p-value           = ', num2str(1-chi2cdf(stat,1)) ]);
    disp(' ');

    % Plot volatility function for alternative values of gam
    tt = seqa(1946+12/12,1/12,t);
    figure(1)

    subplot(2,2,1);
    plot(tt,drt./r1t.^0.0);
    title('$\gamma=0.0$')
    box off
    axis tight

    subplot(2,2,2);
    plot(tt,drt./r1t.^0.5);
    title('$\gamma=0.5$')
    box off
    axis tight
    
    subplot(2,2,3);
    plot(tt,drt./r1t.^1.0);
    title('$\gamma=1.0$')
    box off
    axis tight

    subplot(2,2,4);
    plot(tt,drt./r1t.^1.5);
    title('$\gamma=1.5$')
    box off
    axis tight

end
%
%------------------------- Functions -------------------------------------%
%
%-------------------------------------------------------------------------%
% Define the moment equations 
%-------------------------------------------------------------------------%
function dt = meqn(b,drt,r1t)
    
        ut = drt - b(1) - b(2)*r1t;
        zt = [ones(size(ut,1),1),r1t];
        dt = repmat(ut,1,2).*zt;
        dt = [dt,repmat((ut.^2 - (b(3)^2)*r1t.^(2*b(4)) ),1,2).*zt];
   
end
%-------------------------------------------------------------------------%
% Defines the mean of the moment conditions  
%-------------------------------------------------------------------------%
function ret = meaneqn(b,drt,r1t)

        ret = (mean(meqn(b,drt,r1t)))';

end
%-------------------------------------------------------------------------%
% GMM objective function which also computes the optimal w   
%-------------------------------------------------------------------------%   
function [ret,w] = q(b,drt,r1t,lmax,flag)
        
    t = length(drt);
    d = meqn(b,drt,r1t);
    g = mean(d)';

    if flag
        w   = d'*d;
        tau = 1;
        while tau <= lmax
            wtau = d((tau+1):size(d,1),:)'*d(1:(size(d,1)-tau),:);
            w    = w + (1.0-tau/(lmax+1))*(wtau + wtau');
            tau  = tau + 1;
        end
        w = w./t;            
    else
        w = eye(size(d,2));
    end
        ret = g'*inv(w)*g;
end