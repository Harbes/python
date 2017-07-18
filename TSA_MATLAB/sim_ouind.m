%=========================================================================
%
%   Monte Carlo analysis to investigate the sampling properties 
%   of the indirect estimator of Ornstein-Uhlenbech process.
%
%      Gourieroux et. al. (1993) J of Appl. Eco.
%
%=========================================================================
function sim_ouind( )

    clear all;
    clc;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) );

    % Parameters of Ornstein-Uhlenbech process
    t     = 250;                                    
    kappa = 0.8;
    alpha = 0.1;
    sig2  = 0.06^2;
                 
    y0    = 0.1;            %     Initial value of Ornstein-Uhlenbech process  
    dt    = 0.1;            %     Continuous time step interval       
    h     = 10/dt;          %     Scalar to control the simulation run
    n     = t*h;            %     Length of the simulated series      

    % Simulation settings
    opt   = optimset('LargeScale','off','Display','off');
    nreps = 200;
    theta = zeros(nreps,3);

    % Main DO LOOP to generate sampling disribution
    for j = 1:nreps
    
        disp(['Replication number ', num2str(j) ]);

        % Generate the actual data for the Ornstein-Uhlenbech process       
        u  = randn(t,1);                                                           
        y  = recserar(alpha*(1-exp(-kappa)) + ...
            sqrt(sig2)*sqrt((1-exp(-2*kappa))/(2*kappa))*u,y0,exp(-kappa));
    
        % Estimate the auxiliary model using actual data  
        x  = [ones(size(y,1)-1,1),trimr(y,0,1)];
        y  = trimr(y,1,0);
        b  = x\y;

        bhat    = zeros(3,1);
        bhat(1) = b(1)/(1 - b(2));
        bhat(2) = 1 - b(2);
        bhat(3) = mean((y - x*b).^2);
     
        % Compute the indirect estimator
        v          = sqrt(dt)*randn(n,1);
        b0         = [alpha ; kappa ; sig2];
        b          = fminsearch(@(b) q(b,dt,v,y0,bhat,n),b0,opt);
        b(3)       = abs(b(3));
        theta(j,:) = b';
    end
    
    % Generate statistics on the sampling distribution
    m     = mean(theta);
    stdev = std(theta);
    rmse  = sqrt(mean(bsxfun(@minus,theta,[alpha kappa sig2]).^2));

    disp(' ');
    disp(['Number of replications              =  ', num2str(nreps) ]);
    disp(['Sample size                         =  ', num2str(t) ]);
    
    disp('    True      Mean      Std.err.  RMSE ')
    disp([[alpha kappa sig2]' m' stdev' rmse']);

end
%
%------------------------- Functions -------------------------------------%
%
%-------------------------------------------------------------------------%
% Objective function to compute the indirect estimator 
%-------------------------------------------------------------------------%
function retp = q(b,dt,v,y0,bhat,n)

        b(3) = abs(b(3));                                                    
        ys   = recserar(b(1)*b(2)*dt + sqrt(b(3))*v,y0,1-b(2)*dt);                      
        
        % Generate discrete time data by choosing every t-th observation  
        nys  = zeros(n*dt,1);
        for  i = 1:n*dt
        
            temp   = i*(1/dt); 
            nys(i) = ys(temp);                       
        
        end

        ys = nys;

        xs  = [ones(size(ys,1)-1,1),trimr(ys,0,1)];
        ys  = trimr(ys,1,0);

        b   = xs\ys;

        bhats    = zeros(3,1);
        bhats(1) = b(1)/(1 - b(2));
        bhats(2) = 1 - b(2);
        bhats(3) = mean((ys - xs*b).^2);

        retp = ( (bhat - bhats)'*(bhat - bhats) );

end
        
