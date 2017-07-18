%=========================================================================
%
%   Investigate the effects of misspecifying the trend 
%
%=========================================================================
function nts_trend( )

    clear all
    clc
    
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',15) );

    % DGP: Simulate I(1) variable
    t     = 100;
    delta = 0.5;
    sig2  = 1.0;
    y0    = 0.0;
    y     = recserar(delta + sqrt(sig2)*randn(t,1),y0,1.0);          

    % Estimate the deterministic trend model      
    [b,se] = ar(y,2);       

    disp('    Est      Std Err.   t-stat')
    disp([b se b./se])

	% Sampling distribution of the slope parameter estimate 
    % of the deterministic trend model   
    ndraws = 50000;
    t      = 5000;

    b1 = zeros(ndraws,1);

    for i = 1:ndraws

        % DGP is I(1)
        y      = recserar(delta + sqrt(sig2)*randn(t,1),y0,1.0);                                           
        [b,~] = ar(y,2);                               
        b1(i) = b(2);

    end

    % Compute the standardized statistic   
    stat = sqrt(t)*(b1 - delta);                        

    disp(' ')
    disp(['Mean of scaled statistic     = ', num2str(mean(stat)) ])
    disp(['Mean (theoretical)           = ', num2str(0.0) ])
    disp(['Variance of scaled statistic = ', num2str(std(stat)^2) ])
    disp(['Variance (theoretical)       = ', num2str(sig2*6/5) ])

   
    hist(stat,31);
end

%
%--------------------------- Functions ----------------------------------
% 
%-------------------------------------------------------------------------
% Regression estimates of trend models
%-------------------------------------------------------------------------
function [b,se]=ar(y,ntrend)

    t = length(y);
    % ntrend = 0 (constant and lag)
    if ntrend == 0.0 
        
        x = [ones(t-1,1)   trimr(y,0,1) ];                                  
        y = trimr(y,1,0);  
         
    % ntrend = 1 (constant, time trend and lag)    
    elseif ntrend == 1.0
        x = [ones(t-1,1)  seqa(0,1,t-1)'/100  trimr(y,0,1)];     
        y = trimr(y,1,0);   
    
    % ntrend = 2 (constant and time trend)    
    else
        
        x = [ones(t,1)  seqa(0,1,t)'/1];                          
        y = trimr(y,0,0);   
        
    end    

    b = x\y;                                     
    e = y - x*b;                                 
    s2 = e'*e/length(e);                            
    se = sqrt( diag( s2*inv(x'*x) ) );          

end