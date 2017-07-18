%=========================================================================
%
%   Estimate the ACD model of Engle-Russell for AMR on Aug01, 2006 
%
%=========================================================================
function discrete_acd( )  

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
    
    % Define durations
    y    = data(:,7);
    t    = length(y);                       
    
    % Descriptive statistics
    disp(['Total number of observations             = ',num2str(t) ]);
    disp(['Average duration time (in seconds)       = ',num2str(mean(y)) ]);
    disp(['Minimum duration time (in seconds)       = ',num2str(min(y)) ]);
    disp(['Maximum duration time (in seconds)       = ',num2str(max(y)) ]);
    disp(['Unconditional Hazard rate                = ',num2str(1/mean(y)) ]);
    disp(' ');

    %**********************************************************************
    %***
    %***     Generate graph
    %***
    %********************************************************************** 
    
    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;
    
    %--------------------------------------------------------%
    % Panel (a)
    subplot(1,2,1)
    hist(y,21)
    title('(a) Histogram of duration times $y_t$','FontSize',8);
    box off
    set(gca,'FontSize',8);
    
    
    %--------------------------------------------------------%
    % Panel (b)
    subplot(1,2,2)
    bar(seqa(1,1,50),acf(y,50))
    title('(b) ACF of duration times $y_t$','FontSize',8);
    box off
    set(gca,'FontSize',8);

    % Print the tex file to the relevant directory
    %laprint(1,'durations','options','factory','keepfontprops', 'on');


    % Estimate the ACD model
    theta0 = [  0.9670198086627571 
                0.9323466594671925 
                0.0000091497133080 ];

    theta0 = [0.1 ; 0.1 ; 0.9];
    ops    = optimset('LargeScale','off','Display','off');

    [ thetac,lf ] = fminunc(@(b) neglog(b,y),theta0,ops);
    
    disp(['Log-likelihood (ACD) model    = ',num2str(-lf) ])
    disp( ' ' )

    % Estimate the GARCH model and show the equivalence with the ACD model    

    theta0 = [  0.9670198086627571 
                0.9323466594671925 
                0.0000091497133080 ];

    [ thetag,lf ] = fminunc(@(b) negloggarch(b,y),theta0,ops);

    disp(['Log-likelihood (GARCH) model  = ',num2str(-lf) ])
    disp( ' ' )

    disp(' Parameter Estimates ')
    disp('    ACD       GARCH ')
    disp([ abs(thetac) abs(thetag)] )

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  The log-likelihood of the ACH model 
%-------------------------------------------------------------------------
function lf = neglog(b,y) 

    b = abs(b);
    mu = recserar( (b(1)^2) + (b(2)^2)*trimr( [0.0;y],0,1),mean(y),(b(3)^2));     
    lf = -mean ( -log(mu) - y./mu );

end
%-------------------------------------------------------------------------
%  The log-likelihood of the GARCH model assuming conditional normality 
%-------------------------------------------------------------------------
function lf = negloggarch(b,y)

    u  = sqrt(y);
    s2 = recserar( (b(1)^2) + (b(2)^2)*trimr( [0.0; u.^2],0,1),mean(y),(b(3)^2));     
    z  = u./sqrt(s2);                                                     

    lf = -mean (- 0.5*log(2*pi) - 0.5*log(s2) - 0.5*z.^2 );

end
%--------------------------------------------------------------------------
% ACF function based on running a sequence of OLS regressions  
%--------------------------------------------------------------------------
function acorr = acf(y,lags)

    acorr = zeros(lags,1);
    t     = length(y);

    for i = 1:lags
        
        y0 = trimr(y,i,0);
        x  = [ ones(t-i,1) trimr(y,0,i)];
        b  = inv(x'*x)*x'*y0;
        
        acorr(i) = b(2);
    end
end

