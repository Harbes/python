% ========================================================================
%
%   Monte Carlo analysis to investigate the sampling properties
%   of the indirect estimator of geometric Brownian motion.
%
%   Gourieroux et. al. (1993) J of Appl. Eco.
% ========================================================================
function sim_geobrind( )

    clear all;
    clc;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',12345) );

    t     = 100;            %  Sample size                 
    h     = 1;
    mu    = 1.0;
    sig   = 0.5;
    theta = [mu;sig];
    y0    = 1;
    dt    = 0.01;



    %    Main DO LOOP to generate sampling disribution
    ndraws = 100;
    bmsm   = zeros(ndraws,2);

    for j = 1:ndraws

        disp(['Replication number ', num2str(j) ]);

        % Generate the actual data                        
        yact = exp(recserar(mu - 0.5*sig^2 + sig*randn(t,1),log(y0),1.0));

        % Estimate the auxiliary model using actual data   
        bhat = zeros(2,1);
        bhat(1) = mean(trimr(yact,1,0)./trimr(yact,0,1)) - 1;
        bhat(2) = sqrt(mean((trimr(yact,1,0)./trimr(yact,0,1) - 1 - bhat(1)).^2));
  
        % Generate errors to simulate the true model      
        dw = randn(t/dt,h)*sqrt(dt);

        % Estimate the indirect model using simulated data   
        b0  = [mu;sig];
        opt = optimset('LargeScale','off','Display','off');
        b   = fminunc(@(b) fobj(b,bhat,h,dt,dw,y0,t),b0,opt);
      
        bmsm(j,:) = b';

    end
    % Generate statistics on the sampling distribution    
    m     = mean(bmsm);
    stdev = sqrt(mean((bmsm - repmat(m,size(bmsm,1),1)).^2));
    rmse  = sqrt(mean((bmsm - repmat(theta',size(bmsm,1),1)).^2));

    disp(' ');
    disp(['Number of replications              = ', num2str(ndraws) ]);
    disp(['Sample size                         = ', num2str(t) ]);
    disp(' ');
    disp(['True population parameter           = ', num2str(theta')]);
    disp(['Mean of estimates                   = ', num2str(m) ]);
    disp(['Standard deviation of estimates     = ', num2str(stdev) ]);
    disp(['RMSE of Theta                       = ', num2str(rmse) ]);

end
%
% ------------------------ Functions ------------------------------------%
%
%-------------------------------------------------------------------------
%     The objective function to compute the indirect estimator  
%-------------------------------------------------------------------------
function retp = fobj(b,bhat,h,dt,dw,y0,t)
        
        btilda = zeros(size(bhat,1),h);
        ysim   = exp(recserar( (b(1) - 0.5*b(2)^2)*dt + b(2)*dw,log(y0)*ones(1,h),ones(1,h)) );

        % Generate discrete time data by choosing every t-th observation  
        nys  = zeros(t,1);
        for  i = 1:t
        
            temp   = i*(1/dt); 
            nys(i) = ysim(temp);                       
        
        end

        yt = nys;

        btilda(1,:) = (mean(trimr(yt,1,0)./trimr(yt,0,1)) - 1)';
        btilda(2,:) = sqrt(mean((trimr(yt,1,0)./trimr(yt,0,1) - 1 - btilda(1,:)).^2))';
       
        w    = eye(size(bhat,1));
        retp = (bhat-mean(btilda')')'*w*(bhat - mean(btilda')');

        
end




