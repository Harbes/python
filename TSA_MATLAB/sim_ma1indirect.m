%=========================================================================
%
%      Monte Carlo analysis to investigate the sampling properties
%      of the indirect estimator of a first order MA model.
%
%      Gourieroux et. al. (1993) J of Appl. Eco.
%
%=========================================================================
function sim_ma1indirect( )
    
    clear all;
    clc;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) );

    % Parameters of MA(1) process
    t     = 250;                                
    theta = 0.5;                    
    lag   = 3;             

    % Simulation settings
    opt    = optimset('LargeScale','off','Display','off');
    nreps = 1000;
    b  = zeros(nreps,1);

    % Main DO LOOP to generate sampling disribution
    for j = 1:nreps
    
        % Generate the actual data for the MA(1) process       
        u  = randn(t,1);
        y  = trimr(u,1,0) - theta*trimr(u,0,1);
        
        % Estimate the the auxiliary model using actual data  
        y   = y - mean(y);
        if lag == 1
            bhat = (trimr(y,0,1))\trimr(y,1,0); 
        elseif lag == 2
            bhat = [trimr(y,1,1),trimr(y,0,2)]\trimr(y,2,0); 
        elseif lag == 3 
            bhat = [trimr(y,2,1),trimr(y,1,2),trimr(y,0,3)]\trimr(y,3,0); 
        end
    
        % Compute indirect estimator (could use a line search algorithm)
        e    = randn(t,1);
        b(j) = fminsearch(@(b) q(b,e,lag,bhat),0.5,opt);
        
    end
    
    % Generate statistics on the sampling distribution
    b = tanh(b);
    m     = mean(b);
    stdev = std(b);
    rmse  = sqrt(mean(b-theta)^2);

    disp(' ');
    disp(['Number of replications              =  ', num2str(nreps) ]);
    disp(['Sample size                         =  ', num2str(t) ]);
    
    disp('    True      Mean      Std.err.  RMSE ')
    disp([theta m stdev rmse]);
        
end
%
%------------------------- Functions -------------------------------------%
%
%-------------------------------------------------------------------------%
% Objective function to compute the indirect estimator 
%-------------------------------------------------------------------------%
function val = q(b,e,lag,bhat)

    % Simulate data making sure that b(1) is in the unit circle
    ys     = trimr(e,1,0) - tanh(b)*trimr(e,0,1);        
    ys     = ys - mean(ys);
 
    if lag == 1
        bhats = (trimr(ys,0,1))\trimr(ys,1,0); 
    elseif lag == 2
        bhats = [trimr(ys,1,1),trimr(ys,0,2)]\trimr(ys,2,0); 
    elseif lag == 3
        bhats = [trimr(ys,2,1),trimr(ys,1,2),trimr(ys,0,3)]\trimr(ys,3,0); 
    end
    val = (bhat - bhats)'*(bhat - bhats);      
end
        