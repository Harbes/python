%=========================================================================
%
%      Program to compute the statistical properties of daily equity returns 
%
%=========================================================================

function garch_statistic( )

    clear all
    clc
     
    % Load data  
    load equity

    % Choose equity index and compute percentage return
    equity  = ftse;     % ftse, dow, nikkei                                             
    y       = 100*(trimr( log( equity ),1,0 ) - trimr( log( equity ),0,1 )); 
       
    % Compute the autocorrelation of returns and returns squared
    mlag  = 20;
    acfy  = zeros( mlag+1,1 );
    acfy2 = zeros( mlag+1,1 );
    
    acfy(1)  = acfunction( y,0 );
    acfy2(1) = acfunction( y.^2,0 );
    
    for j = 2:mlag 
    
        acfy(j)  = acfunction( y,j );
        acfy2(j) = acfunction( y.^2,j );
    end
  
    % Compute the moments of the unconditional (empirical) distribution
    y_z = (y - mean(y))/std(y);
    disp( ' Skewness ');
    disp( mean(y_z.^3) );
    disp( ' Kurtosis ');
    disp( mean(y_z.^4) );

    %*********************************************************************
    %**     Generate graphs
    %*********************************************************************
    
    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;

    tt   = 1:1:length(y);
    lags = 0:1:mlag;
    
    subplot(2,2,1)
    plot(tt,y,'-k','LineWidth',0.75);
    title('(a) Returns');
    ylabel('$y_t$');
    xlabel('$t$');
    box off;
    axis tight;

    subplot(2,2,2)
    plot(tt,y.^2,'-k','LineWidth',0.75);
    title('(b) Squared Returns');
    ylabel('$y_t^2$');
    xlabel('$t$');
    box off;
    axis tight;


    subplot(2,2,3)
    %plot( lags,acfy,'-k','LineWidth',0.75);
    bar( lags,acfy,0.65,'k' )
    title('(c) ACF Returns');
    ylabel('$\text{acf}(y_t)$');
    xlabel('$p$');
    box off;
    axis tight;

    subplot(2,2,4)
    %plot( lags,acfy2,'-k','LineWidth',0.75);
    bar( lags,acfy2,0.65,'k' )
    title('(d) ACF Squared Returns');
    ylabel('$\text{acf}(y_t^2)$');
    xlabel('$p$');
    box off;
    axis tight;
   
    %laprint(1,'fig-volatility1','options','factory');

    % Kernel density graph
    figure(2);
    clf;
    
    xi = -6:0.1:6;
    f  = ksdensity( y,xi);
    ft = normpdf( xi,0,1 );
    
    plot( xi,f','-k',xi,ft','--k','LineWidth',0.75 );
    ylabel('$f(y)$');
    xlabel('$y$');
    
    %laprint(2,'fig-volatility2','options','factory');
end
%
%--------------------------- Functions ----------------------------------
% 
%-------------------------------------------------------------------------
%   Compute the ACF at kag p
%-------------------------------------------------------------------------

function r = acfunction( x,p )

    d = x - mean( x );
    r = sum( trimr(d,p,0).*trimr(d,0,p) ) / sum(d.^2);

end




