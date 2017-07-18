%=========================================================================
%
%   Sample moments of stochastic and deterministic trend models
%
%=========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',15) );

%  Choose model
% itrend = 0 for stochastic trend model 
% itrend = 1 for deterministic trend model
itrend = 0;

% Parameters
t      = 800;
ndraws = 50000;
sig2   = 1.0;
y0     = 0.0;

% Parameters of the stochastic trend model
delta = 0.0;                
phi   = 1.0;     % 0.8 for stationary case Example 16.4

%  Parameters of the deterministic trend model
beta1 = 0.2;
beta0 = 0.1;                

% Initialise arrays
m1 = zeros(ndraws,1);       
m2 = zeros(ndraws,1);        
m3 = zeros(ndraws,1);        
m4 = zeros(ndraws,1);       

for i = 1:ndraws

    if itrend == 0                          

        % AR(1) model with first observation discarded
        y    = trimr( recserar(delta + sqrt(sig2)*randn(t+1,1),y0,phi) , 1 , 0);           
        ylag = trimr(y,0,1);

        if phi < 1.0
            
            % Standardization based on t for stationary case
            m1(i) = sum(ylag.^1)/t^1;                                                            
            m2(i) = sum(ylag.^2)/t^1;                    
            m3(i) = sum(ylag.^3)/t^1;                
            m4(i) = sum(ylag.^4)/t^1;
        else
            
            % Standardization for nonstationary case
            m1(i) = sum(ylag.^1)/t^1.5;                                                              
            m2(i) = sum(ylag.^2)/t^2;                       
            m3(i) = sum(ylag.^3)/t^2.5;                  
            m4(i) = sum(ylag.^4)/t^3;
        end

    elseif itrend == 1               

        trend = seqa(1,1,t);
        y     = beta0 + beta1*trend + sqrt(sig2)*randn(t,1);                

        % Standardization based on t for stationary case
        m1(i) = sum(trend.^1)/t^2;                                                               
        m2(i) = sum(trend.^2)/t^3;
        m3(i) = sum(trend.^3)/t^4;
        m4(i) = sum(trend.^4)/t^5;

    end

end


disp(['Sample size    = ',num2str(t) ])
disp(' ')
disp(['Mean of m1     = ',num2str(mean(m1)) ])
disp(['Variance of m1 = ',num2str(std(m1)^2) ])
disp(' ')
disp(['Mean of m2     = ',num2str(mean(m2)) ])
disp(['Variance of m2 = ',num2str(std(m2)^2) ])
disp(' ')
disp(['Mean of m3     = ',num2str(mean(m3)) ])
disp(['Variance of m3 = ',num2str(std(m3)^2) ])
disp(' ')
disp(['Mean of m4     = ',num2str(mean(m4)) ])
disp(['Variance of m4 = ',num2str(std(m4)^2) ])


