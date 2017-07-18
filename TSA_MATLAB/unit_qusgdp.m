%=========================================================================
%
%   Program to illustrate basic concepts using U.S. GDP
%
%=========================================================================
function unit_qusgdp( )

    clear all
    clc

    % Load the data
    load usgdp.mat

    t  = length(gdp);
    y = log(gdp);
    
        
    %**********************************************************************
    %***     Generate graph
    %**********************************************************************

    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;
  
    dvec = seqa(1947,1/4,t);

    %--------------------------------------------------------%
    plot(dvec,y,'-k');
    box off
    ylim( [7.0 10 ] );
    set(gca,'ytick',[7.0 7.5 8.0 8.5 9.0 9.5 10.0]);
    xlim( [ 1940 2010 ] );
    set(gca,'xtick',[1950 1960 1970 1980 1990 2000 2010]);
    ylabel('Log of Real GDP');
    xlabel('Year')
    
    % Print the tex file to the relevant directory
    %laprint(1,'usgpd','options','factory');
      
    % Detrending
    x = [ ones(t,1) seqa(1,1,t)' ];
    
    [ uols,bols ] = glsdetrend(y,x,-t);
    [ udif,bdif ] = glsdetrend(y,x,0);
    [ ugls,bgls ] = glsdetrend(y,x,-13.5);
        
    disp('Constant and trend - OLS detrending')
    disp(bols')
    disp('Constant and trend - Differencing')
    disp(bdif')
    disp('Constant and trend - GLS detrending')
    disp(bgls')
    disp('')
    
    %**********************************************************************
    %***     Generate graph
    %**********************************************************************
    figure(2);
    clf;

    %--------------------------------------------------------%
    plot(dvec,uols,'-k', ...
         dvec,udif,'--k',...
         dvec,ugls,':k');
    box off
    ylim( [-0.15 0.15 ] );
    set(gca,'ytick',[-0.15 -0.10 -0.05 0.00 0.05 0.10 0.15]);
    xlim([ 1940 2010 ] );
    set(gca,'xtick',[1950 1970 1990 2010]);

    ylabel('Residuals');
    xlabel('Year')

    
    % Print the tex file to the relevant directory
    %laprint(2,'usgpdresids','options','factory');

    
    % Perform unit root tests on GDP 
    % Assuming no autocorrelation
    k = 0;
    
    % Assuming autocorrelation
    %k = kmaic(uols);
    disp(['Number of lags to adjust for autocorrelation = ',num2str(k) ])
        
    % adf tests based on ols detrending
    adf_ols = adftests(uols,k);
    disp(' ');
    disp('OLS Detrending')
    disp('  DFa         DFt     ');
    disp(adf_ols');
    
    % mtests based on ols detrending 
    mtests_ols = mtests(uols,k);
    disp(' ');
    disp('OLS Detrending')
    disp('  mza         msb       mzt     ');
    disp(mtests_ols);

    % adf tests based on gls detrending
    adf_gls = adftests(ugls,k);
    disp(' ');
    disp('OLS Detrending')
    disp('  DFa         DFt     ');
    disp(adf_gls');
    
    % mtests based on gls detrending 
    mtests_gls = mtests(ugls,k);
    disp(' ');
    disp('GLS Detrending')               
    disp('  mza         msb       mzt     ');
    disp(mtests_gls);

    % OLS with trend break 1973Q1
    dt   = cumsum( dvec > 1973 )';

    xd   = [ x  dt ]; 
    bd   = xd\y;          
    yhat = xd*bd;
    ud   = y - yhat;
 
    % Linear trend  no break at 1973Q1
    yhat1 = xd*([bd(1:2) ; 0] );

    %**********************************************************************
    %***     Generate graph
    %**********************************************************************
    figure(3);
    clf;

    %--------------------------------------------------------%
    subplot(1,2,1)
    plot(dvec,y,'-k', ...
         dvec,yhat,'--k',...
         dvec,yhat1,':k');
 
    title('(a)')
    ylim( [7.0 10 ] );
    set(gca,'ytick',[7.0 7.5 8.0 8.5 9.0 9.5 10.0]);
    xlim( [ 1940 2010 ])
    set(gca,'xtick',[1950 1970 1990 2010]);
    ylabel('Time Trends, Real GDP');
    xlabel('Year')

    %--------------------------------------------------------%
    subplot(1,2,2)
    plot(dvec,uols,'-k', ...
         dvec,ud,'--k');
    title('b')
    ylim( [-0.15 0.15 ] );
    set(gca,'ytick',[-0.15 -0.10 -0.05 0.00 0.05 0.10 0.15]);
    xlim( [ 1940 2010 ])
    set(gca,'xtick',[1950 1970 1990 2010]);
    ylabel('Residuals');
    xlabel('Year')

    % Print the tex file to the relevant directory
    %laprint(3,'usgpdtrends','options','factory');

    % Unit root tests with break at 1973q1
    uols   = ud;
    k      = kmaic(uols);
    ugls   = glsdetrend(y,xd,-18.4);
    MZtols = trimr(mtests(uols,k)',2,0);
    MZtgls = trimr(mtests(ugls,k)',2,0);

    disp('Mzt tests based on break in 1973q1 ');
    disp(['OLS Detrending    ', num2str(MZtols) ])                
    disp(['GLS Detrending    ', num2str(MZtgls) ])              

    %  Estimate break point from levels OLS
    TL = floor(0.15*t);
    ssel = std(y)^2*(t-1)*ones(t,1);

    for j = TL:t-TL 
  
        x = [ones(t,1) seqa(1,1,t)' [ zeros(j,1) ; seqa(1,1,t-j)'] ];
        e = y-x*(x\y);
        
        ssel(j) = e'*e/length(e);
    end
    
    [~,TBl] = min(ssel);
        
    x    = [ones(t,1) seqa(1,1,t)' [zeros(TBl,1);seqa(1,1,t-TBl)']];
    ytbl = x*(x\y);

    % Re-estimate break point from GLS using cbar from levels OLS breakpoint
    taul = TBl/t;
    cbar = -18.5;
    sseqd = std(y)^2*(t-1)*ones(t,1);
    
    for j = TL:t-TL 
  
        x  = [ones(t,1) seqa(1,1,t)' [zeros(j,1); seqa(1,1,t-j)']];
        yc = [ y(1); trimr(y,1,0)-(1+cbar/length(y))*trimr(y,0,1)];
        xc = [x(1,:) ; trimr(x,1,0)-(1+cbar/length(y))*trimr(x,0,1)];
        e  = yc-xc*(xc\yc);
  
        sseqd(j) = e'*e/length(e);
    end
    [~,TBqd] = min(sseqd);

    x     = [ones(t,1) seqa(1,1,t)' [zeros(TBqd,1) ; seqa(1,1,t-TBqd)']];
    ytbqd = x*(x\y);
 
    % Unit root tests using GLS breakpoint
    cbar   = -18.4;
    uols   = glsdetrend(y,x,-t);
    k      = kmaic(uols);
    ugls   = glsdetrend(y,x,cbar);
    MZtols = trimr(mtests(uols,k)',2,0);
    MZtgls = trimr(mtests(ugls,k)',2,0);

    disp('Mzt tests based on GLS estimate of break ');
    disp(['OLS Detrending    ', num2str(MZtols) ])                
    disp(['GLS Detrending    ', num2str(MZtgls) ])              

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
%  M tests
%-------------------------------------------------------------------------
function tests = mtests(u,k)

    s2  = ar1rvar(u,k);
    n   = length(u);
    tmp = sum(u(1:n-1).^2);
    
    %disp(['Sum of u(1) to u(t-1) = ',num2str(tmp) ])
    %disp(['Last value: u(t)      = ',num2str(u(n)) ])
    %disp(' ') 

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

