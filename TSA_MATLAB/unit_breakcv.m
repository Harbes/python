%=========================================================================
%
%   5% cvs for unit root tests with break in trend at known point
%
%=========================================================================
function unit_breakcv( )
 
    clear all
    clc
    
    reps = 1000;

    Tv   = [25,50,100,250,500,1000];
    tauv = seqa(0.15,0.05,15);
    cbar = [-17.6,-17.8,-18.2,-18.4,-18.6,-18.4,-18.4,-18.2,...
            -18.0,-17.6,-17.4,-17.0,-16.6,-16.0,-15.2];

    % Allocate arrays
    dfols    = zeros(reps,length(tauv)); 
    dfgls    = dfols; 
    mztols   = dfols;
    dfolscv  = zeros(length(Tv),length(tauv)); 
    dfglscv  = dfolscv; 
    mztolscv = dfolscv;
   
    for Tc = 1:length(Tv)

        t = Tv(Tc);
        disp(['Current sample size = ',num2str(t) ])

        for tauc = 1:length(tauv)
      
            TB = floor(tauv(tauc)*t);
            DT = [zeros(TB,1); seqa(1,1,t-TB)'];
            x  = [ones(t,1) seqa(1,1,t)' DT];
    
            RandStream.setDefaultStream( RandStream('mt19937ar','seed',42) );
        
            for rep=1:reps 
                y = cumsum(randn(t,1));
                uols = glsdetrend(y,x,-t);
                ugls = glsdetrend(y,x,cbar(Tc));
                dfols(rep,tauc) = trimr(adftests(uols,0),1,0);
                dfgls(rep,tauc) = trimr(adftests(ugls,0),1,0);
                mztols(rep,tauc) = trimr(mtests(uols,0)',2,0);
            end
        end
        dfolscv(Tc,:) = quantile(dfols,0.05);
        dfglscv(Tc,:) = quantile(dfgls,0.05);
        mztolscv(Tc,:) = quantile(mztols,0.05);

    end
    format short g
    
    disp('DF-OLS')
    disp('    tauv')
    disp([ Tv' dfolscv ])

    disp('')
    disp('DF-GLS')
    disp('    tauv')
    disp([ Tv' dfglscv ])
    
    disp('')
    disp('MZt-OLS')
    disp('    tauv')
    disp([ Tv' mztolscv ])

end
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
