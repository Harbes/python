% ========================================================================
%
%   Estimate a multifactor model of the Stock-Watson business cycle model
%   using indirect estimation based on Gallant and Tauchen.
%   Auxiliary model is a VAR.
%
% ========================================================================
function sim_bcycle( )

    clear all
    clc
    
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',12345) );


    % Parameter values
    t = 600;                                
    lam = [ 1 ; 0.8 ; 0.7 ; 0.5 ; -0.5 ; -1 ];  
    sig = 0.2*seqa(1,1,6)';
    phi = 0.8;                              
    
    % Generate data
    eta = randn(t,1);
    z   = randn(t,6);
    s   = recserar( eta, 0.0, phi);
    y   = bsxfun(@times,lam',s) + bsxfun(@times,sig',z);


    % Estimate the auxiliary model by a sequence of ar(1) regressions  
    g = [ ];
    bhat = [ ];

    for i = 1:6;
        yvar = trimr(y(:,i),2,0);
        xvar = [ones(length(yvar),1)  trimr(y(:,i),1,1)  trimr(y(:,i),0,2) ];
        b    = xvar\yvar;
        vhat = yvar - xvar*b;
        s2   = mean(vhat.^2);
        bhat = [ bhat  [b ; s2] ];
        
        % Moment conditions based on the data
        g = [ g  [bsxfun(@times,vhat,xvar)/s2   (0.5*vhat.^2/s2^2 - 0.5/s2)] ];    

    end
    iinv = inv(g'*g/length(g));           

    % Simulation Estimation  
    n    = 10*t;                  % Length of simulation run 
    etas = randn(n,1);            % Fix simulated disturbances   
    zs   = randn(n,6);

    % Call optimiser
    opt    = optimset('LargeScale','off','Display','off');
    theta0 = [ lam'  sig'  phi] ;

	[theta,qmin] = fminunc(@(b) q(b,bhat,etas,zs,iinv),theta0,opt);

    dof  = length(bhat(:))-length(theta);
    jstat = t*qmin;
    
    disp('Parameter estimates');
    disp('  True   Estimated ')
    disp([theta0' theta'])
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
    phi = b(13);
        
    ys = bsxfun(@times,lam,recserar( etas, 0.0, phi)) + bsxfun(@times,sig,zs);
        
    gs = [ ];
    for i = 1:6
        
        % Estimate auxiliary model evaluated at bhat  
        yvar  = trimr(ys(:,i),2,0);          
        xvar  = [ones(length(yvar),1) trimr(ys(:,i),1,1) trimr(ys(:,i),0,2)];
        vhats = yvar - xvar*bhat([1 2 3],i);
        gs    = [gs  bsxfun(@times,vhats,xvar)/bhat(4,i)  (0.5*vhats.^2/bhat(4,i)^2 - 0.5/bhat(4,i)) ];
    end

    ret = mean(gs)*iinv*mean(gs)';   

end