%=========================================================================
%
%    Program to test the "Great Moderation" hypothesis, using 
%    real annual US GDP data per capita from 1946 to 2006.
%                                                         
%=========================================================================

function hetero_moderation( ) 	


    clear all;
    clc;

    
    gdp = [ 11241;
            10924;
            11206;
            10957;
            11717;
            12412;
            12668;
            13032;
            12719;
            13389;
            13410;
            13435;
            13088;
            13783;
            13840;
            13933;
            14552;
            14971;
            15624;
            16420;
            17290;
            17532;
            18196;
            18573;
            18392;
            18771;
            19555;
            20485;
            20195;
            19961;
            20822;
            21565;
            22526;
            22982;
            22666;
            23007;
            22347;
            23146;
            24593;
            25382;
            26024;
            26664;
            27514;
            28221;
            28429;
            28007;
            28556;
            28941;
            29741;
            30128;
            30880;
            31886;
            32833;
            33904;
            34755;
            34645;
            34837;
            35361;
            36300;
            37052;
            37752 ];

   
    % Compute growth rate
    y = 100*( log(gdp(2:end)) - log(gdp(1:end-1)) );
  
%**************************************************************************
%**
%**     Generate graph
%**
%**************************************************************************
    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;

    plot(1947:1:2006,y,'-k')
    xlabel('Years');
    ylabel('Growth Rate of US Per Capital GDP');
    axis tight;
    box off;

    %     Generate descriptive statistics   

    disp(['Sample mean (1947-1983)     = ', num2str(mean(y(1:37))) ]);
    disp(['Sample mean (1984-2006)     = ', num2str(mean(y(38:60))) ]);
    disp(['Sample variance (1984-2006) = ', num2str(mean((y(1:37)-mean(y(1:37))).^2)) ]);
    disp(['Sample variance (1984-2006) = ', num2str(mean((y(38:60)-mean(y(38:60))).^2)) ]);

    d = [zeros(37,1) ; ones(23,1)];       %     Construct dummy variable   
    t = length(y);
   
    % Estimate the unconstrained model 
    theta = 0.1*ones(4,1);
    flag  = 1;
    [theta1,a1,a,a,a,h1] = fminunc(@(b) neglog(b,y,d,flag),theta);  
   
    omega1 = inv(h1);
    lnl1 = -a1;                         

    disp(['ML estimate of the variance (1946-1983) = ',num2str(exp(theta1(3))) ]);
    disp(['ML estimate of the variance (1984-2006) = ',num2str(exp(theta1(3) + theta1(4))) ]);

    % Estimate the constrained model    
    theta = 0.1*ones(3,1);
    [aaa,a0] = fminunc(@(b) neglog0(b,y,d,flag),theta);        
    lnl0=-a0;
 
    % LR test for heteroskedasticity
    lr = -2*t*(lnl0 - lnl1);
    disp(['LR statistic            = ',num2str(lr) ]); 
    disp(['p-value                 = ',num2str(1-cdf('chi2',lr,1)) ]);

    % Wald test for heteroskedasticity
    r = [0 , 0 , 0 , 1];
    q = 0;
    wd = t*(r*theta1 - q)'*inv(r*omega1*r')*(r*theta1 - q);
    disp(['Wald statistic          = ',num2str(wd) ]); 
    disp(['p-value                 = ',num2str(1-cdf('chi2',wd,1)) ]);

    % LM test (regression form)    
    x = [ones(t,1) , d];   
    % Stage 1 regression
    b = x\y; 
    u = y - x*b;    
    w = [ones(t,1) , d];                            
    v = u.^2;
    % Stage 2 regression
    b  = w\v; 
    e  = v - w*b;
    r2 = 1 - sum(e.^2)/sum( (v-mean(v)).^2 );
    lm = t*r2;
    disp(['LM statistic            = ',num2str(lm) ]); 
    disp(['p-value                 = ',num2str(1-cdf('chi2',lm,1)) ]);

    % Wald test that beta1 = 0  
    r  = [0 , 1 , 0 , 0];
    q  = 0;
    wd = (r*theta1 - q)'*inv(r*omega1*r')*(r*theta1 - q);
 
    disp(' ')
    disp(['Wald statistic of beta1 = 0   = ',num2str(wd) ]);
    disp(['p-value                       = ',num2str(1 - cdf('chi2',wd,1)) ]);

  %-----------------------------------------------------------------------
  % Tests based on the assumption that beta1 = 0.
  %-----------------------------------------------------------------------
  
    % Estimate the unconstrained model but with beta1 = 0     
    theta = [0.1 ; 0.1 ; 0.1];
    flag = 0;
    [theta1,a1,a,a,a,h1] = fminunc(@(b) neglog(b,y,d,flag),theta);  
   
    omega1 = inv(h1);
    lnl1 = -a1;      
   
    % Estimate the constrained model but with beta1 = 0   
    theta = [0.1 ; 0.1 ];
    [aa,a0] = fminunc(@(b) neglog0(b,y,d,flag),theta);        
    
    omega1 = inv(h1);
    lnl0=-a0;
   
    % LR test for heteroskedasticity
    lr = -2*t*(lnl0 - lnl1);
    disp(['LR statistic            = ',num2str(lr) ]); 
    disp(['p-value                 = ',num2str(1-cdf('chi2',lr,1)) ]);

    % Wald test   
    r = [0 , 0 , 1];
    q = 0;
    wd = t*(r*theta1 - q)'*inv(r*omega1*r')*(r*theta1 - q);
    disp(['Wald statistic          = ',num2str(wd) ]); 
    disp(['p-value                 = ',num2str(1-cdf('chi2',wd,1)) ]);
 
    % LM test (regression form)   
    x = ones(t,1);                                
    % Stage 1 regression
    b = x\y; 
    u = y - x*b;    
    w = [ones(t,1) , d];                            
    v = u.^2;
    % Stage 2 regression
    b  = w\v; 
    e  = v - w*b;
    r2 = 1 - sum(e.^2)/sum( (v-mean(v)).^2 );
    lm = t*r2;
    disp(['LM statistic            = ',num2str(lm) ]); 
    disp(['p-value                 = ',num2str(1-cdf('chi2',lm,1)) ]);

end


%
%-------------------------Functions------------------------------------
%
%-----------------------------------------------------------------------
%      Negative log-likelihood function (unconstrained)   
%-----------------------------------------------------------------------
function lf = neglog( b,y,x,flag )

    if flag
        mu   = b(1) + b(2)*x;
        sig2 = exp(b(3) + b(4)*x);
    else
        mu = b(1);
        sig2 = exp(b(2) + b(3)*x);
    end
    
    lnl  = -(1/2)*log(2*pi*sig2) - (y - mu).^2 ./(2*sig2); 
    lf = -mean( lnl );
    
end
%-----------------------------------------------------------------------
%      Negative log-likelihood function (constrained)   
%-----------------------------------------------------------------------
function lf = neglog0( b,y,x,flag )

    if flag
        mu   = b(1) + b(2)*x;
        sig2 = exp(b(3) + 0*x);
    else
        mu   = b(1);
        sig2 = exp(b(2) + 0*x);
    end
       
    lnl  = -(1/2)*log(2*pi*sig2) - (y - mu).^2 ./(2*sig2);
    lf   = -mean( lnl );
    
end


%     Define the log of the likelihood for model with beta1 = 0 (unconstrained)         
 
function lnl_t = lnltbo(b,d,y)

 
     muet = b(1) + 0*d;
     sigt = sqrt(exp(b(3) + b(4)*d));
     lnl_t = -(1/2)*log(2*pi*sigt.^2) - (y - muet).^2 ./(2*sigt.^2);        % Log like at each obs    
end

%     Define the constrained log of the likelihood for model with beta1 = 0
%     
function lnl_t = lnltbo0(b,d,y);

 
     muet = b(1) + 0*d;
     sigt = sqrt(exp(b(3) + 0*d));
     lnl_t = -(1/2)*log(2*pi*sigt.^2) - (y - muet).^2 ./(2*sigt.^2);        % Log like at each obs    
end
