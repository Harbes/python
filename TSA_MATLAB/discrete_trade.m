%=========================================================================
%
%   Estimate the ACH model of Hamilton-Jorda for AMR on Aug01, 2006 
%
%   This program can't reproduce the GAUSS results reported in the book
%   because the likelihood function is flat relative to 
%   the delta1 parameter in particular. Consequently the results
%   are particularly sensitive to choice of algorithm and starting values.
%
%=========================================================================
function discrete_trade( )  

    clear all
    clc

    % Load data for AMR on Aug 1, 2006.
    % Order of the variables:
    %   1.  Hour
    %   2.	Minute
    %   3.	Second
    %   4.	Time (time effectively starts at 9.30am and ends at 4.00pm with the time interval being one second)
    %   5.	y (a binary variable giving a one if a trade has occured, and 0 otherwise
    %   6.	N (a counter variable which increases by one when a trade occurs)
    %   7.	u (gives the duration between last trade)
    %   8.	Firm News dummy (AMR)
    %   9.	Macro News dummy

    load amr_aug1.mat;

    % Trim the first 32 observations and exclude the last observation
    % Ensure postive durations
    data = trimr(data,32,1);
    
    % Define variables
    y     = data(:,5);
    nt    = data(:,6);
    u     = data(:,7);
    firm  = data(:,8);
    macro = data(:,9);
    
    ubar      = mean(u);           % Average duration time                      
    firm_lag  = trimr(firm,0,1);   % Lagged firm news variable                              
    macro_lag = trimr(macro,0,1);  % Lagged macroeconomic news variable                         
    y         = trimr(y,1,0);      % Lineup other data  
    u         = trimr(u,1,0);
    t         = length(y);

    % Descriptive statistics
    disp(['Total number of observations             = ',num2str(t) ]);
    disp(['Total number of trades                   = ',num2str(sum(y)) ]);
    disp(['Average number of trades (per second)    = ',num2str(mean(y)) ]);
    disp(['Average duration time (in seconds)       = ',num2str(mean(u)) ]);
    disp(['Minimum duration time (in seconds)       = ',num2str(min(u)) ]);
    disp(['Maximum duration time (in seconds)       = ',num2str(max(u)) ]);
    disp(['Unconditional Hazard rate                = ',num2str(1/mean(u)) ]);
    disp(' ');


    % Estimate the restricted model (without news variables)
    ops = optimset('LargeScale','off','Display','off');
    theta_0 =   [   0.0196493517339459 
                    0.0010464918437577 
                    0.9840278397571428 ];
    [ theta0,l0 ] = fminunc(@(b) neglog0(b,y,u,ubar),theta_0,ops);
    
    disp('Restricted parameter estimates')
    disp(theta0)
    
    l0 = -l0;   
    disp(['Log-likelihood (restricted) = ',num2str(l0) ]);
    disp(['TxLog-likelihood        = ',num2str(t*l0) ]);
    disp(' ');


    % Compute h at optimat parameters
    [ h,si ] = hazard0(theta0,u,ubar); 
    disp('ACH Model without annoucement variables');
    disp(['Hazard rate                         = ',num2str(mean(h)) ]);
    disp(['Time to the next trade (in seconds) = ',num2str(1/mean(h)) ]);
    disp(' ');

	% Estimate the unrestricted model (with news variables)  
    ops = optimset('LargeScale','off','Display','off');
    theta_0 =   [     theta0(1)
                      theta0(2)
                      theta0(3)
                      0.1
                      0.1               ];
    
   [ theta1,l1 ] = fminsearch(@(b) neglog(b,y,u,ubar,firm_lag,macro_lag),theta_0,ops);
   
    disp('Unrestricted parameters')
    disp(theta1)
    l1 = -l1;   
    disp(['Log-likelihood (unrestricted) = ',num2str(l1) ]);
    disp(['TxLog-likelihood              = ',num2str(t*l1) ]);
    disp(' ');

    % Compute h at optimal parameters
    [ h,si ] = hazard(theta1,u,ubar,firm_lag,macro_lag); 
    disp('ACH Model without annoucement variables');
    disp(['Hazard rate                         = ',num2str(mean(h)) ]);
    disp(['Time to the next trade (in seconds) = ',num2str(1/mean(h)) ]);
    disp(' ');


    disp('ACH Model without annoucement variables: general model');
    h_adj = 1/(1 + exp(mean(si)));         
    disp(['Hazard rate                         = ',num2str(h_adj) ]);
    disp(['Time to the next trade (in seconds) = ',num2str(1/h_adj) ]);
    disp(' ');

    disp('ACH Model without firm news variables');
    h_adj = 1./(1 + exp(mean(si) + theta1(5)));    
    disp(['Hazard rate                         = ',num2str(h_adj) ]);
    disp(['Time to the next trade (in seconds) = ',num2str(1/h_adj) ]);
    disp(' ');


    disp('ACH Model without macro news variables');
    h_adj = 1/(1 + exp(mean(si) + theta1(4)));    
    disp(['Hazard rate                         = ',num2str(h_adj) ]);
    disp(['Time to the next trade (in seconds) = ',num2str(1/h_adj) ]);
    disp(' ');

    disp('ACH Model with both firm and macro news variables');
    h_adj = 1./(1 + exp(mean(si) + theta1(4) + theta1(5)));    
    disp(['Hazard rate                         = ',num2str(h_adj) ]);
    disp(['Time to the next trade (in seconds) = ',num2str(1/h_adj) ]);
    disp(' ');

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  The log-likelihood of the restricted ACH model
%-------------------------------------------------------------------------
function lf = neglog0(b,y,u,ubar)

     % ACD(1,1) model
     si = recserar( b(1) + b(2)*trimr([0.0; u],0,1),ubar,b(3));     

     % Hazard rate   
     h = 1./(1 + exp(si));      

     lf = - mean( y.*log(h) + (1 - y).*log(1 - h) );

end

%-------------------------------------------------------------------------
%  Return h and psi for the restricted ACH model
%-------------------------------------------------------------------------
function [h,si] = hazard0(b,u,ubar)

     % ACD(1,1) model
     si = recserar( b(1) + b(2)*trimr([0.0; u],0,1),ubar,b(3));     

     % Hazard rate   
     h = 1./(1 + exp(si));      

end

%-------------------------------------------------------------------------
%  The log-likelihood of the unrestricted ACH model
%-------------------------------------------------------------------------
function lf = neglog(b,y,u,ubar,firm_lag,macro_lag)

     % ACD(1,1) model
     si = recserar( b(1) + b(2)*trimr([0.0; u],0,1),ubar,b(3));     

     % Hazard rate   
     h = 1./(1 + exp(si + b(4)*firm_lag + b(5)*macro_lag));      

     lf = -mean(  y.*log(h) + (1 - y).*log(1 - h) );

end
%-------------------------------------------------------------------------
%  Return h and si for the unrestricted ACH model
%-------------------------------------------------------------------------
function [ h,si ]  = hazard(b,u,ubar,firm_lag,macro_lag)

     % ACD(1,1) model
     si = recserar( b(1) + b(2)*trimr([0.0; u],0,1),ubar,b(3));     

     % Hazard rate   
     h = 1./(1 + exp(si + b(4)*firm_lag + b(5)*macro_lag));      

end

