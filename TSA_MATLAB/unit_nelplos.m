%=========================================================================
%
%   Unit root tests: Nelson-Plosser data (1860 to 1970).
%
%=========================================================================          
function unit_nelplos

    clear all;
    clc;

    % Load Nelson-Plosser data set
    load nelson_plosser.mat

    % Variable 'data' constains the following variables
    % Date 1860 - 1970
    % Real GNP, 1909 to 1970  
    % Nominal GNP, 1909 to 1970
    % Real per capita GNP, 1909 to 1970
    % Industrial production, 1860 to 1970
    % Employment, 1890 to 1970
    % Unemployment rate, 1890 to 1970
    % GNP deflator, 1889 to 1970
    % Consumer prices index, 1860 to 1970
    % Wages, 1900 to 1970
    % Real wages, 1900 to 1970 
    % Money stock, 1889 to 1970
    % Velocity, 1869 to 1970
    % Bond yield, 1900 to 1970
    % SP500, 1871 to 1970

    % Take logs except for bond yield
    rgnp   = log(data(50:111,2));
    gnp    = log(data(50:111,3));
    pcrgnp = log(data(50:111,4));
    ip     = log(data(1:111,5));
    emp    = log(data(31:111,6));
    un     = log(data(31:111,7));
    prgnp  = log(data(30:111,8));
    cpi    = log(data(1:111,9));
    wg     = log(data(41:111,10));
    rwg    = log(data(41:111,11));
    m      = log(data(30:111,12));
    vel    = log(data(10:111,13));
    bnd    = data(41:111,14);
    sp500  = log(data(12:111,15));
    
    % Choose data and lag length and flag unemployment)
    y    = un;
    flag = 1;
    lags = 6;
    
    
    [acy,acdy,acu,Mztols,Mztgls] = np(y,lags,flag);
    
    disp('         Autocorrelation Functions')
    disp('      Lag     Levels    1st Diff. Detrend')
    disp('    -------------------------------------')
    disp( [seqa(1,1,lags)' acy acdy acu] )
    
    disp('     Unit Root Tests ')
    disp('    Mztols     Mztgls')
    disp('    ------------------')
    disp( [Mztols Mztgls] )
    
    % Union of rejections
    if ~flag
        
        cvgls = -2.867; 
        cvols = -3.130; 
        tau = 1.038;
    else
        cvgls = -2.867; 
        cvols = -3.130; 
        tau = 1.038;
    end
        
    disp(['UOR Mztols classification is I(',num2str(Mztols > tau*cvols),')'])
    disp(['UOR Mztgls classification is I(',num2str(Mztgls > tau*cvols),')'])

end
%--------------------------- Functions -----------------------------------
% 
%--------------------------------------------------------------------------
%   Compute acf and Mtests for Nelson Plosser data  
%--------------------------------------------------------------------------
function [acy,acdy,acu,Mztols,Mztgls] = np(y,lags,flag)

    t = length(y);
    
    if ~flag
        % acf of y
        acy = acf(y,lags);
    
        % acf acf of differences
        acdy = acf(trimr(y,1,0)-trimr(y,0,1),lags); 
    
        % acf of residuals from a linear trend
        x    = [ones(length(y),1) seqa(1,1,t)']; 
        uols = glsdetrend(y,x,-t);
        acu  = acf(uols,lags);

        % MZt tests
        cbar   =-13.5; 
        k      = kmaic(uols);
        ugls   = glsdetrend(y,x,cbar);
        Mztols = trimr(Mtests(uols,k)',2,0);
        Mztgls = trimr(Mtests(ugls,k)',2,0);
    
    % Unemployment flagged
    else
        % acf of y
        acy = acf(y,lags);
    
        % acf acf of differences
        acdy = acf(trimr(y,1,0)-trimr(y,0,1),lags); 
    
        % acf of residuals from a constant only
        x    = ones(length(y),1); 
        uols = glsdetrend(y,x,-t);
        acu  = acf(uols,lags);

        % MZt tests
        cbar   =-7; 
        k      = kmaic(uols);
        ugls   = glsdetrend(y,x,cbar);
        Mztols = trimr(Mtests(uols,k)',2,0);
        Mztgls = trimr(Mtests(ugls,k)',2,0);       
    end
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
function tests = Mtests(u,k)

    s2  = ar1rvar(u,k);
    n   = length(u);
    tmp = sum(u(1:n-1).^2);
    
    disp(['Sum of u(1) to u(t-1) = ',num2str(tmp) ])
    disp(['Last value: u(t)      = ',num2str(u(n)) ])
    disp(' ') 

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
