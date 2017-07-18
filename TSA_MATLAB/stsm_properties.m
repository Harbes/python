%=========================================================================
%
%   Program to identify properties of US macro variables
%
%=========================================================================
function stsm_properties( )

    clear all
    clc

    format short
    % Read the data: quarterly US data from Jan-1959 to Dec-1998
    load sims_data.mat
 
    % Define varaibles
    r    = ytdata(:,1);        
    lex  = log( ytdata(:,2) );
    lcp  = log( ytdata(:,3) );
    lm   = log( ytdata(:,4) );
    lp   = log( ytdata(:,5) );
    lo   = log( ytdata(:,6) );
    sdum = ytdata(:,7:17);
    
    % Number of lags used to compute the ACF and PACF
    lags = 6; 
    
    % Compute the ACF and PACF on levels of variables
    y     = [ r lm lp lo ];
    acf1  = acf( y(:,1),lags );
    pacf1 = pacf( y(:,1),lags );
    acf2  = acf( y(:,2),lags );
    pacf2 = pacf( y(:,2),lags );
    acf3  = acf( y(:,3),lags );
    pacf3 = pacf( y(:,3),lags );
    acf4  = acf( y(:,4),lags );
    pacf4 = pacf( y(:,4),lags );
   
    disp(' ACF of series in levels ')
    disp('-------------------------')
    disp('      Lag       y1       y2        y3        y4')
    disp( [ (1:1:lags)' acf1 acf2 acf3 acf4 ] );
    
    disp(' PACF of series in levels ')
    disp('-------------------------')
    disp('      Lag       y1       y2        y3        y4')
    disp( [ (1:1:lags)' pacf1 pacf2 pacf3 pacf4 ] );
    
    
    % Compute the ACF and PACF on first difference of variables
    dy   = trimr(y,1,0) - trimr(y,0,1);         
    acf1  = acf( dy(:,1),lags );
    pacf1 = pacf( dy(:,1),lags );
    acf2  = acf( dy(:,2),lags );
    pacf2 = pacf( dy(:,2),lags );
    acf3  = acf( dy(:,3),lags );
    pacf3 = pacf( dy(:,3),lags );
    acf4  = acf( dy(:,4),lags );
    pacf4 = pacf( dy(:,4),lags );
   
    disp(' ACF of first difference of the series ')
    disp('-----------------------------------------')
    disp('      Lag       y1       y2        y3        y4')
    disp( [ (1:1:lags)' acf1 acf2 acf3 acf4 ] );
    
    disp(' PACF of first difference of the series ')
    disp('-----------------------------------------')
    disp('      Lag       y1       y2        y3        y4')
    disp( [ (1:1:lags)' pacf1 pacf2 pacf3 pacf4 ] );
    
    % Compute the ACF and PACF on second difference of variables
    ddy = trimr(y,12,0) - trimr(y,0,12);       
    acf1  = acf( ddy(:,1),lags );
    pacf1 = pacf( ddy(:,1),lags );
    acf2  = acf( ddy(:,2),lags );
    pacf2 = pacf( ddy(:,2),lags );
    acf3  = acf( ddy(:,3),lags );
    pacf3 = pacf( ddy(:,3),lags );
    acf4  = acf( ddy(:,4),lags );
    pacf4 = pacf( ddy(:,4),lags );
   
    disp(' ACF of 12th difference of the series ')
    disp('-----------------------------------------')
    disp('      Lag       y1       y2        y3        y4')
    disp( [ (1:1:lags)' acf1 acf2 acf3 acf4 ] );
    
    disp(' PACF of 12th difference of the series ')
    disp('-----------------------------------------')
    disp('      Lag       y1       y2        y3        y4')
    disp( [ (1:1:lags)' pacf1 pacf2 pacf3 pacf4 ] );
end

%
%--------------------------- Functions  ----------------------------------
% 
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
%--------------------------------------------------------------------------
% PACF function based on running a sequence of OLS regressions  
%--------------------------------------------------------------------------
function pcorr = pacf(y,lags)

    pcorr = zeros(lags,1);
    t     = length(y);
    x     = ones(t,1);

    for i = 1:lags
        
        y0 = trimr(y,i,0);
        x = [trimr(x,1,0) trimr(y,0,i)];
        b = inv(x'*x)*x'*y0;
        
        pcorr(i) = b(i+1);
    end
end


