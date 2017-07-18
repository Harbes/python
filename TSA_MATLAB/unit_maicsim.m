%=========================================================================
%
%   Finite sample properties of tests and autocorrelation
%
%=========================================================================
function unit_maicsim( ) 

    clear all
    clc
    
    t     = 200;
    reps  = 10000;
    trend = 1;
    
    if trend
        x     = [ones(t,1) seqa(1,1,t)']; 
        cbar  = -13.5; 
        MZtcv = -2.87;
    else
        x = ones(t,1); 
        cbar = -7; 
        MZtcv = -1.94;
    end
    
    arv = [0,0.3,0.6]; 
    %arv = 0; 
    %mav = 0;
    mav = [-0.8,-0.4,0,0.4,0.8];
    %phiv = {1,0.9,0.8,0.7,0.6,0.5};
    phiv = [1,0.95,0.9,0.85];

    Maic  = zeros(length(phiv),length(arv),length(mav));
    Mmaic = zeros(length(phiv),length(arv),length(mav));
    Maic  = zeros(length(phiv),length(arv),length(mav));
    DFgls = zeros(length(phiv),length(arv),length(mav));

    for phic = 1:length(phiv) 
        for arc = 1:length(arv)
            for mac =1:length(mav) 
                for rep = 1:reps
          
                    u = randn(t+1,1);
                    u = trimr(u,1,0)-mav(mac)*trimr(u,0,1);
                    u = recserar(u,u(1),arv(arc));
                    y = recserar(u,u(1),phiv(phic));
            
                    % Detrending
                    uols = glsdetrend(y,x,-t); 
                    ugls = glsdetrend(y,x,cbar); 
            
                    % Select lag length using MAIC
                    k = kmaic(uols); 
        
                    Mmaic(phic,arc,mac) = Mmaic(phic,arc,mac) + (trimr(mtests(ugls,k)',2,0)<MZtcv)/reps;
                    Maic(phic,arc,mac)  = Maic(phic,arc,mac) + (trimr(mtests(ugls,kaic(uols))',2,0)<MZtcv)/reps;
                    DFgls(phic,arc,mac) = DFgls(phic,arc,mac)+ (trimr(adftests(ugls,k),1,0)<MZtcv)/reps;
                
                end
                disp('')
                disp(['phi          =   ',num2str(phiv(phic)) ])
                disp(['rho          =   ',num2str(arv(arc))  ])
                disp(['theta        =   ',num2str(mav(mac))  ])
                disp('-----------------------------------')
                disp(['MZt-GLS-maic  =   ',num2str(Mmaic(phic,arc,mac)) ])
                disp(['MZt-GLS-aic  =    ',num2str(Maic(phic,arc,mac)) ])
                disp(['DFt-GLS-maic  =   ',num2str(DFgls(phic,arc,mac)) ])
                disp('')
                
            end
        end
    end
end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Detrending function: 
%       cbar = -7 constant; 
%       cbar = -13.5 linear trend  
%       cbar = -T for OLS detrending
%-------------------------------------------------------------------------
function [ u,b ] = glsdetrend( y,x,cbar )

    t = length( y );
    yc = [ y(1)  ; (trimr(y,1,0)-(1+cbar/t)*trimr(y,0,1)) ];
    xc = [x(1,:) ; (trimr(x,1,0)-(1+cbar/t)*trimr(x,0,1)) ];
    
    b  = xc\yc;
    u  = y - x*b;
 
 
end
%-------------------------------------------------------------------------
%  ADF coefficient and t tests: u must be already de-trended
%-------------------------------------------------------------------------
function adfstats= adftests(u,k)


    du = trimr(u,1,0) - trimr(u,0,1);
    x  = trimr(u,0,1);

    % Set up lags of du
    if k>0;
    
        ldu = lagmatrix(du,1:1:k);
        x   = trimr( [x ldu],k,0 );
    end
    
    xxi = inv(x'*x); 
    b   = xxi*x'*trimr(du,k,0); 
    e   = trimr(du,k,0)-x*b; 
    s2  = e'*e/length(e);

    adfstats = [length(u)*b(1) ; b(1)/sqrt(s2*xxi(1,1))];

end
%-------------------------------------------------------------------------
%  Select ADF lag length by AIC: u should be OLS de-trended residuals
%-------------------------------------------------------------------------
function k = kaic(u)

    kmax = floor(12*(length(u)/100)^.25); 
    aic = zeros(kmax+1,1);

    % Set up lagged regressors
    du = trimr(u,1,0)-trimr(u,0,1);
    x  = [ trimr(u,0,1) lagmatrix(du,seqa(1,1,kmax)) ];
    x(any(isnan(x),2),:) = [];

    for j =0:kmax
  
        b = x(:,1:1+j)\trimr(du,kmax,0); 
        e = trimr(du,kmax,0)-x(:,1:1+j)*b; 
        s2 = e'*e/length(e);
  
        aic(j+1) = log(s2) + 2*(j+1)/length(e);
    end
    
    [~,k] = min(aic);
    k     = k-1;
    
end

%-------------------------------------------------------------------------
%  Select ADF lag length by MAIC: u should be OLS de-trended residuals
%-------------------------------------------------------------------------
function k = kmaic(u)

    kmax = floor(12*(length(u)/100)^.25); 
    maic = zeros(kmax+1,1);

    % Set up lagged regressors
    du = trimr(u,1,0)-trimr(u,0,1);
    x  = [ trimr(u,0,1) lagmatrix(du,seqa(1,1,kmax)) ];
    x(any(isnan(x),2),:) = [];

    for j =0:kmax
  
        b = x(:,1:1+j)\trimr(du,kmax,0); 
        e = trimr(du,kmax,0)-x(:,1:1+j)*b; 
        s2 = e'*e/length(e);
  
        maic(j+1) = log(s2) + 2*(j+b(1)^2*sum(x(:,1).^2)/s2)/length(e);
    end
    
    [~,k] = min(maic);
    k     = k-1;
    
end
%-------------------------------------------------------------------------
%  M tests
%-------------------------------------------------------------------------
function tests = mtests(u,k)

    s2  = ar1rvar(u,k);
    n   = length(u);
    tmp = sum(u(1:n-1).^2);
    
    u2  = tmp/n^2;     
    mza = (u(n)^2/n-s2)/(2*u2);     
    msb = sqrt(u2/s2);
    mzt = msb*mza;
    
    tests = [mza msb mzt];
 
end
%-------------------------------------------------------------------------
%  Autoregressive long run variance estimator
%-------------------------------------------------------------------------
function s2 = ar1rvar(u,k)

    du = trimr(u,1,0)-trimr(u,0,1); 
    x  = trimr(u,0,1);

     if k > 0 
        x  = [ x lagmatrix(du,seqa(1,1,k)) ];
        x(any(isnan(x),2),:) = [];
     end


    b = x\trimr(du,k,0); 
    e = trimr(du,k,0) - x*b; 
    s2 = e'*e/length(e);

    if k>0 
        s2 = s2/(1-sum(trimr(b,1,0)))^2; 
    end

end

