% ========================================================================
%
%   Estimate a multifactor model of the Stock-Watson business cycle model
%   using indirect estimation based on Gallant and Tauchen.
%   Auxiliary model is a VAR.
%
% ========================================================================
function sim_stockwatson( )

    clear all
    clc
    
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',12345) );

    % Load monthly data September 1959 to September 2009 for Australia
    load('bcycle.mat')
    employed = data(:,1);    % Employed     
    gdp      = data(:,2);    % GDP          
    hincome  = data(:,3);    % Household income   
    indprod  = data(:,4);    % Industrial production     
    retail   = data(:,5);    % Retail sales  
    unemp    = data(:,6);    % Unemployment rate
    index    = data(:,7);    % Coincident index 

    y = log( [employed gdp hincome indprod retail 1./unemp] );
    y = 100*(trimr(y,3,0) - trimr(y,0,3));
    y = bsxfun(@minus,y,mean(y));
    y = bsxfun(@rdivide,y,std(y));
    
    t = length(y);
    
    % corrx((trimr(ln(index),3,0)-trimr(ln(index),0,3))~y);
    % Estimate the auxiliary model by a sequence of ar(1) regressions  
    g = [ ];
    bhat = [ ];

    for i = 1:6;
        yvar = trimr(y(:,i),2,0);
        xvar = [trimr(y(:,i),1,1)  trimr(y(:,i),0,2) ];
        b    = xvar\yvar;
        uhat = yvar - xvar*b;
        s2   = mean(uhat.^2);
        bhat = [ bhat  [b ; s2] ];
        
        % Moment conditions based on the data
        g = [ g  [bsxfun(@times,uhat,xvar)/s2   (0.5*uhat.^2/s2^2 - 0.5/s2)] ];    
 
    end
    iinv = inv(g'*g/length(g));           

    % Simulation Estimation  
    n    = 30*t;                  % Length of simulation run 
    etas = randn(n,1);            % Fix simulated disturbances   
    zs   = randn(n,6);

    % Call optimiser
    opt          = optimset('LargeScale','off','Display','iter');           
    theta0       = [rand(1,12) 1];
	[theta,qmin] = fminunc(@(b) q(b,bhat,etas,zs,iinv),theta0,opt);
    theta(end)   = tanh(theta(end));
    dof  = length(bhat(:))-length(theta);
    jstat = t*qmin;
    
    disp('Parameter estimates');
    disp('Parameter Estimates ')
    disp(theta')
    disp(['Value of the objective function (Q) = ' num2str(qmin) ]);
    disp(['Value of the J-test (TQ)            = ' num2str(jstat) ]);
    disp(['Degrees of freedom                  = ' num2str(dof) ]);
    disp(['P-value                             = ' num2str(1-chi2cdf(jstat,dof)) ]);

end

%
%------------------------- Functions -------------------------------------%
%
%-------------------------------------------------------------------------%
% Objective function to compute the EMM estimator 
%-------------------------------------------------------------------------%
function  ret = q(b,bhat,etas,zs,iinv)

    lam = b(1:6);
    sig = b(7:12);
    phi = tanh(b(13));
        
    ys = bsxfun(@times,lam,recserar( etas, 0.0, phi)) + bsxfun(@times,sig,zs);
    ys = bsxfun(@minus,ys,mean(ys));
    
    gs = [ ];
    for i = 1:6
        
        % Estimate auxiliary model evaluated at bhat  
        yvar  = trimr(ys(:,i),2,0);          
        xvar  = [ trimr(ys(:,i),1,1) trimr(ys(:,i),0,2)];
        vhats = yvar - xvar*bhat([1 2],i);
        gs    = [gs  bsxfun(@times,vhats,xvar)/bhat(3,i)  (0.5*vhats.^2/bhat(3,i)^2 - 0.5/bhat(3,i)) ];
    end

    ret = mean(gs)*iinv*mean(gs)';   

end
