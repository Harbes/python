% ========================================================================
%
%      This program performs a Monte Carlo analysis to investigate the
%      sampling properties of the Indirect estimator of Brownian motion with drift.
%
%      Gourieroux et. al. (1993) J of Appl. Eco.
% ========================================================================

function sim_brownind( )
    clear all;
    clc;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',12345) );

    t     = 500;            %     Sample size                         
    mu    = 0.5;
    sig2  = 0.25;
    theta = [mu;sig2];      %     Parameters of Brownian motion       
    y0    = 1;              %     Initial value of Brownian process   
    dt    = 0.1;            %     Continuous time step interval       
    h     = 1/dt;           %     Scalar to control the simulation run
    n     = t*h;            %     Length of the simulated series      


    %    Main DO LOOP to generate sampling disribution
    ndraws     = 100;
    theta_ind  = zeros(ndraws,2);

    for j = 1:ndraws

        disp(['Replication number ', num2str(j) ]);

        %  Simulate sample data using an exact discretisation      
        u  = randn(t,1);                                                     
        y  = recserar( mu + sqrt(sig2)*u , y0 , 1.0 );                            

        % Estimate the auxiliary model using actual data    
        dy   = trimr(y,1,0) - trimr(y,0,1);                
        bhat = [mean(dy) ; mean((dy - mean(dy)).^2)];

        % Generate errors to be used to compute the emm estimator.
        % Note that these errors need to be generated outside of the procedure.
        v = sqrt(dt)*randn(n,1);

        % Estimate model  
        bstart  = [mu ; sig2];
        options = optimset('LargeScale','off','Display','off'); 
        b       = fminunc(@(b) q(b,dt,v,y0,bhat,n),bstart,options);
        
        theta_ind(j,:) = b';

    end


    % Generate statistics on the sampling distribution    
    mean_ = mean(theta_ind);
    stdev = sqrt(mean((theta_ind - repmat(mean_,size(theta_ind,1),1)).^2));
    rmse  = sqrt(mean((theta_ind - repmat(theta',size(theta_ind,1),1)).^2));

    disp(' ');
    disp(['Number of replications              = ', num2str(ndraws) ]);
    disp(['Sample size                         = ', num2str(t) ]);
    disp(' ');
    disp(['True population parameter           = ', num2str(theta')]);
    disp(['Mean of estimates                   = ', num2str(mean_) ]);
    disp(['Standard deviation of estimates     = ', num2str(stdev) ]);
    disp(['RMSE of Theta                       = ', num2str(rmse) ]);

end
%
% ------------------------ Functions ------------------------------------%
%
%-------------------------------------------------------------------------
%     The objective function to compute the indirect estimator  
%-------------------------------------------------------------------------
function val = q(b,dt,v,y0,bhat,n)

        b(2)= abs(b(2));                                
        
        % Generate continuous time data and choose every t-th observation  
        ys  = recserar( b(1)*dt + sqrt(b(2))*v , y0 , 1 );                                               
        nys = zeros(n*dt,1);
       
        for  i = 1:n*dt
        
            temp   = i*(1/dt); 
            nys(i) = ys(temp);                       
        
        end      
        ys    = nys;
        dys   = trimr(ys,1,0) - trimr(ys,0,1);
        bhats = [mean(dys) ; mean((dys - mean(dys)).^2)];
        val   = (bhat - bhats)'*(bhat - bhats) ;

end



