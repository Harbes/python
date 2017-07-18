%=========================================================================
%
%   A heteroskedastic model of money shocks in US asset markets
%
%=========================================================================
function hetero_event( )

    clear all
    clc
    
    % Load data
    load usasset.mat

    %   1.	tcm3m           (from FRB webpage)
	%	2.	tcm6m  
	%	3.	tcm1y  
	%	4.	tcm2y 
	%	5.	tcm3y
	%	6.	tcm5y
	%	7.	tcm7y
    %   8.  tcm10y
    %   9.  edm1m
    %   10. xrate ($US/$AU)
	%	11.	frb target in US
    %   12. change in target variable US
    %   13. dummy variable corresponding to event dates US (taken from Poole and Kuttner)
    %   14. value-weighted returns in US    (from CRSP)
    %   15. value-weighted index US         (from CRSP)

    % Different to Poole dating
    ydata(1274,13) = 1.0;
    ydata(1275,13) = 0.0;

    % Define variables
    yield    = ydata(:,[1 2 3 4 5 6 7 8]);  % 3m, 6m 1y, 2y, 3y, 5y, 7y, 10y as percentages      
    euro     = ydata(:,9);                  % 1m  euro rate as  percentage                        
    event    = ydata(:,13);                 % dummy variable on meeting dates                                 
    nyse     = ydata(:,15);                 % NYSE equity price index                                         

    % Dependent variables
    dyield = (trimr(yield,1,0) - trimr(yield,0,1));            
    deuro  = trimr(euro,1,0) - trimr(euro,0,1);                  
    rnyse  = 100*(trimr(log(nyse),1,0) - trimr(log(nyse),0,1));   

    y = dyield(:,1);      
    d = trimr(event,1,0);   % Event day
    x = [ones(length(y),1)  d];
    t = length(y);

    % Compute descriptive statistics on non-event and event days    
    y_nonevent = y(d == 0);
    y_event    = y(d == 1);
    
    disp(['Sample mean (event)         = ',num2str(mean(y_event)) ]);
    disp(['Sample mean (non-event)     = ',num2str(mean(y_nonevent)) ]);
    disp(['Sample variance (event)     = ',num2str(mean((y_event-mean(y_event)).^2)) ]);
    disp(['Sample variance (non-event) = ',num2str(mean((y_nonevent-mean(y_nonevent)).^2)) ]);

   %  Estimate the unrestricted model by MLE with starting values based on OLS
    bols = x\y;
    u    = y - x*bols;
    s2   = mean(u.^2); 

    start = [ bols ; s2 ; 0 ];
    ops   = optimset('LargeScale','off','Display','iter');
 
    [bhat1,lf1,~,~,~,hess] = fminunc(@(b) neglog1(b,y,d),start,ops);

    lf1 = -lf1;
    vc1 = (1/t)*inv(hess);
    disp(' ');
    disp('Unrestricted model')
    disp(['MLE estimate of the mean       (event days)     = ',num2str((bhat1(1) + bhat1(2))) ]);
    disp(['MLE estimate of the mean       (non-event days) = ',num2str(bhat1(1)) ]);
    disp(['MLE estimate of the volatility (event days)     = ',num2str(exp(bhat1(3) + bhat1(4))) ]);
    disp(['MLE estimate of the volatility (non-event days) = ',num2str(exp(bhat1(3))) ]);
    disp(' ');

   %  Estimate the restricted model by MLE with starting values based on OLS
    start = [ bols ; s2 ];
    ops   = optimset('LargeScale','off','Display','iter');
 
    [bhat0,lf0,~,~,~,hess] = fminunc(@(b) neglog0(b,y,d),start,ops);
    
    lf0 = -lf0;
    vc0 = (1/t)*inv(hess);

    disp(' ');
    disp('Restricted model')
    disp(['MLE estimate of the mean       (event days)     = ',num2str((bhat0(1) + bhat0(2))) ]);
    disp(['MLE estimate of the mean       (non-event days) = ',num2str(bhat0(1)) ]);
    disp(['MLE estimate of the volatility (event days)     = ',num2str(exp(bhat0(3))) ]);
    disp(['MLE estimate of the volatility (non-event days) = ',num2str(exp(bhat0(3))) ]);
    disp(' ');

    % LR test   
    lr = -2*t*(lf0 - lf1);
    disp(['LR statistic            = ',num2str(lr) ]); 
    disp(['p-value                 = ',num2str(1-cdf('chi2',lr,1)) ]);

    % Wald test 
    r = [0 , 0 , 0 , 1];
    q = 0;
    wd = (r*bhat1 - q)'*inv(r*vc1*r')*(r*bhat1 - q);
    disp(['Wald statistic          = ',num2str(wd) ]); 
    disp(['p-value                 = ',num2str(1-cdf('chi2',wd,1)) ]);

    % LM test (regression form)                                   
    % Stage 1 regression
    b = x\y; 
    u = y - x*b;    
    w = x;                               
    v = u.^2;
    % Stage 2 regression
    b  = w\v; 
    e  = v - w*b;
    r2 = 1 - sum(e.^2)/sum( (v-mean(v)).^2 );
    lm = t*r2;
    disp(['LM statistic (regression) = ',num2str(lm) ]); 
    disp(['p-value                   = ',num2str(1-cdf('chi2',lm,1)) ]);


    
end
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
% Unrestricted log-likelihood function
%-------------------------------------------------------------------------
function lf = neglog1(b,y,d)              

    m  = b(1) + b(2)*d;         % Mean
    s2 = exp(b(3) + b(4)*d);    % Variance
    
    lt = -0.5*log(2*pi) -0.5*log(s2) - 0.5*((y - m).^2)./s2;        
    lf = -mean(lt);

end

%-------------------------------------------------------------------------
% Restricted log-likelihood function
%-------------------------------------------------------------------------
function lf = neglog0(b,y,d)              

    m  = b(1) + b(2)*d;         % Mean
    s2 = exp(b(3));             % Variance
    
    lt = -0.5*log(2*pi) -0.5*log(s2) - 0.5*((y - m).^2)./s2;        
    lf = -mean(lt);

end

    
% 
% 
% 
% 
% /**     Define the log of the likelihood at each observation            **/
% 
% 
% 
% 
