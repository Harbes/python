%=========================================================================
%
%      Over-identification test of the MA(1) model using the EMM estimates
%
%=========================================================================
function sim_ma1overid( )
    
    clear all;
    clc;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',12345) );

 
    % Parameters of MA(1) process
    t     = 250;        % Sample size                         
    theta = 0.5;        % MA1 Population parameter            
    lag   = 3;          % Choose AR lag for auxiliary model   
    p     = 0;          % Used to construct weighting matrix

    % Simulation settings
    opt    = optimset('LargeScale','off','Display','off');
        
    % Generate the actual data for the MA(1) process       
    u  = randn(t,1);
    y  = trimr(u,1,0) - theta*trimr(u,0,1);
        
    % Estimate the the auxiliary model using actual data  
    y   = y - mean(y);
    if lag == 1
        x    = trimr(y,0,1);
        bhat = x\trimr(y,1,0);
        ehat = trimr(y,1,0) - x*bhat;
        g    = ehat.*x;
    elseif lag == 2
        x    = [trimr(y,1,1) , trimr(y,0,2)];
        bhat = x\trimr(y,2,0);
        ehat = trimr(y,2,0) - x*bhat;
        g    = repmat(ehat,1,2).*x;
    elseif lag == 3 
        x    = [trimr(y,2,1) , trimr(y,1,2) , trimr(y,0,3)];
        bhat = x\trimr(y,3,0);
        ehat = trimr(y,3,0) - x*bhat;
        g    = repmat(ehat,1,3).*x;
    end
        
    % Compute the optimal weighting matrix
    i = g'*g;
    
    if p > 0
            
        for l = 1:p
            gam = g((l+1):size(g,1),:)'*g(1:(size(g,1)-l),:);		
            i   = i + (1.0 - l/(p+1))*(gam + gam');			    
        end
    end
    i    = i/length(g);									
    iinv = inv(i);		

    % Compute EMM estimator (could use a line search algorithm)
    e        = randn(t,1);
    [b,qmin] = fminsearch(@(b) q(b,e,lag,bhat,iinv),0.5,opt);
    
    % Number of restrictions  
    r = length(bhat) - length(b);       

    disp(' ');
    disp(['Sample size                         =  ', num2str(t) ]);
    disp(['True parameter value                =  ', num2str(theta) ]);
    disp(['Lags in auxiliary model             =  ', num2str(lag) ]);
    disp(['Over-identification test            =  ', num2str(t*qmin) ]);
    disp(['Number of restrictions              =  ', num2str(r) ]);
    if r > 0
        disp(['p value                             =  ', num2str(1-cdf('chi2',t*qmin,r)) ]);
    end
end
%
%------------------------- Functions -------------------------------------%
%
%-------------------------------------------------------------------------%
% Objective function to compute the EMM estimator 
%-------------------------------------------------------------------------%
 function retp = q(b,e,lag,bhat,iinv)

    ys     = trimr(e,1,0) - tanh(b)*trimr(e,0,1);        
    ys     = ys - mean(ys);
    
    if lag == 1
        xs    = trimr(ys,0,1);
        ehats = trimr(ys,1,0) - xs*bhat;
        gs    = mean( ehats.*xs )';
    elseif lag == 2
        xs    = [trimr(ys,1,1) , trimr(ys,0,2)];
        ehats = trimr(ys,2,0) - xs*bhat;
        gs    = mean( repmat(ehats,1,2).*xs )';
    elseif lag == 3
        xs    = [trimr(ys,2,1) , trimr(ys,1,2) , trimr(ys,0,3)];
        ehats = trimr(ys,3,0) - xs*bhat;
        gs    = mean( repmat(ehats,1,3).*xs )';
    end           
    retp = gs'*iinv*gs ;
 end