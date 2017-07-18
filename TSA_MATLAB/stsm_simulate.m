%=========================================================================
%
%   Simulate an ARMA(2,2) model and compute the ACF and the PACF
%
%==========================================================================
function stsm_simulate( )

    clear all
    clc
    
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',12) );       
    
    t    = 200;                 % Define the sample size      
    lags = 10;                  % Number of lags in the ACF and the PACF      

    % Generate error process

    ut = sqrt(0.1)*randn(t+101,1);      %An additional 100 observations 

    % ARMA(0,0) 
    mue  = 0.0;  
    phi1 = 0.0;    
    phi2 = 0.0;     
    si1  = 0.0;      
    si2  = 0.0;   

    yt    = recserar( mue + trimr(ut,2,0) + si1*trimr(ut,1,1) ...
            + si2*trimr(ut,0,2) , [ 0.0; 0.0 ] , [ phi1;phi2 ] );   
    y0t   = trimr(yt,100,0);                               
    acf0  = acf( y0t,lags );
    pacf0 = pacf( y0t,lags );

   
    % ARMA(1,0)
    mue  = 0.0;  
    phi1 = 0.7;    
    phi2 = 0.0;     
    si1  = 0.0;      
    si2  = 0.0;   

    yt    = recserar( mue + trimr(ut,2,0) + si1*trimr(ut,1,1) ...
            + si2*trimr(ut,0,2) , [ 0.0; 0.0 ] , [ phi1;phi2 ] );   
    y1t   = trimr(yt,100,0);            % Trim the data  
    acf1  = acf( y1t,lags );            % Compute the ACF                             
    pacf1 = pacf( y1t,lags );           % Compute the PACF                            

    
    % ARMA(2,0)
    mue  = 0.0;  
    phi1 = 0.7;    
    phi2 = -0.5;     
    si1  = 0.0;      
    si2  = 0.0;   

    yt    = recserar( mue + trimr(ut,2,0) + si1*trimr(ut,1,1) ...
            + si2*trimr(ut,0,2) , [ 0.0; 0.0 ] , [ phi1;phi2 ] );   
    y2t   = trimr(yt,100,0);          % Trim the data  
    acf2  = acf( y2t,lags );            % Compute the ACF                             
    pacf2 = pacf( y2t,lags );           % Compute the PACF                            

    
    
    % ARMA(0,1)
    mue  = 0.0;  
    phi1 = 0.0;    
    phi2 = 0.0;     
    si1  = 0.9;      
    si2  = 0.0;   

    yt    = recserar( mue + trimr(ut,2,0) + si1*trimr(ut,1,1) ...
            + si2*trimr(ut,0,2) , [ 0.0; 0.0 ] , [ phi1;phi2 ] );   
    y3t   = trimr(yt,100,0);          	
    acf3  = acf(y3t,lags);            % Compute the ACF
    pacf3 = pacf(y3t,lags);           % Compute the PACF 


    % ARMA(0,2)		
    mue  = 0.0;  
    phi1 = 0.0;    
    phi2 = 0.0;     
    si1  = -0.2;      
    si2  = 0.7;   

    yt    = recserar( mue + trimr(ut,2,0) + si1*trimr(ut,1,1) ...
            + si2*trimr(ut,0,2) , [ 0.0; 0.0 ] , [ phi1;phi2 ] );   
    y4t   = trimr(yt,100,0);         
    acf4  = acf(y4t,lags);            % Compute the ACF
    pacf4 = pacf(y4t,lags);           % Compute the PACF 

    % ARMA(1,1)		
    mue  = 0.0;  
    phi1 = 0.8;    
    phi2 = 0.0;     
    si1  = 0.7;      
    si2  = 0.0;   

    yt    = recserar( mue + trimr(ut,2,0) + si1*trimr(ut,1,1) ...
            + si2*trimr(ut,0,2) , [ 0.0; 0.0 ] , [ phi1;phi2 ] );   
    y5t   = trimr(yt,100,0);         
    acf5  = acf(y5t,lags);            % Compute the ACF
    pacf5 = pacf(y5t,lags);           % Compute the PACF 

    
    

    %**********************************************************************
    %***
    %***     Plot the series
    %***
    %**********************************************************************


    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;
    
    t= 0:lags;
     
   
    %--------------------------------------------------------%
    % Panel (a)
    subplot(3,3,1);
    plot(y2t,'-k');
    title('(a) ARMA(2,0)');
    %ylabel('$y_t$');
    %xlabel('$t$');
    axis tight
    box off;

    
    %--------------------------------------------------------%
    % Panel (b)
    subplot(3,3,2);
    bar(t,[1; acf2]);
    title('(b) ACF ARMA(2,0)');
    %ylabel('ACF');
    %xlabel('lag');
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    axis tight
    box off;
     

    %--------------------------------------------------------%
    % Panel (c)
    subplot(3,3,3);
    bar(t,[1; pacf2]);
    title('(b) PACF ARMA(2,0)');
    %ylabel('PACF');
    %xlabel('lag');
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    axis tight
    box off;
 
    %--------------------------------------------------------%
    %--------------------------------------------------------%
    % Panel (d)
    subplot(3,3,4);
    plot(y4t,'-k');
    title('(a) ARMA(0,2)');
    %ylabel('$y_t$');
    %xlabel('$t$');
    axis tight
    box off;


    %--------------------------------------------------------%
    % Panel (e)
    subplot(3,3,5);
    bar(t,[1; acf4]);
    title('(b) ACF ARMA(0,2)');
    %ylabel('ACF');
    %xlabel('lag');
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    axis tight
    box off;
     

    %--------------------------------------------------------%
    % Panel (f)
    subplot(3,3,6);    
    bar(t,[1; pacf4]);
    title('(b) PACF ARMA(0,2)');
    %ylabel('PACF');
    %xlabel('lag');
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    axis tight
    box off;

    %--------------------------------------------------------%
    %--------------------------------------------------------%
    % Panel (g)
    subplot(3,3,7);
    plot(y5t,'-k');
    title('(a) ARMA(1,1)');
    %ylabel('$y_t$');
    xlabel('$t$');
    axis tight
    box off;


    %--------------------------------------------------------%
    % Panel (h)
    subplot(3,3,8);
    bar(t,[1; acf5]);
    title('(b) ACF ARMA(1,1)');
    %ylabel('ACF');
    xlabel('lag');
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    axis tight
    box off;
     

    %--------------------------------------------------------%
    % Panel (i)
    subplot(3,3,9);    
    bar(t,[1; pacf5]);
    title('(b) PACF ARMA(1,1)');
    %ylabel('PACF');
    xlabel('lag');
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    axis tight
    box off;

    
    laprint(1,'stsmsimulate','options','factory');
    
    figure(2);
    clf;
    

    %--------------------------------------------------------%
    % Panel (a)
    subplot(3,3,1);
    plot(y0t,'-k');
    title('(a) ARMA(0,0)');
    %ylabel('$y_t$');
    xlabel('$t$');
    axis tight
    box off;

    
    %--------------------------------------------------------%
    % Panel (b)
    subplot(3,3,2);
    bar(t,[1; acf0]);
    title('(b) ACF ARMA(0,0)');
    %ylabel('ACF');
    xlabel('lag');
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    axis tight
    box off;
     

    %--------------------------------------------------------%
    % Panel (c)
    subplot(3,3,3);
    bar(t,[1; pacf0]);
    title('(b) PACF ARMA(0,0)');
    %ylabel('PACF');
    xlabel('lag');
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    axis tight
    box off;
 
    %--------------------------------------------------------%
    %--------------------------------------------------------%
    % Panel (d)
    subplot(3,3,4);
    plot(y1t,'-k');
    title('(a) ARMA(1,0)');
    %ylabel('$y_t$');
    xlabel('$t$');
    axis tight
    box off;


    %--------------------------------------------------------%
    % Panel (e)
    subplot(3,3,5);
    bar(t,[1; acf1]);
    title('(b) ACF ARMA(1,0)');
    %ylabel('ACF');
    xlabel('lag');
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    axis tight
    box off;
     

    %--------------------------------------------------------%
    % Panel (f)
    subplot(3,3,6);    
    bar(t,[1; pacf1]);
    title('(b) PACF ARMA(1,0)');
    %ylabel('PACF');
    xlabel('lag');
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    axis tight
    box off;

    %--------------------------------------------------------%
    %--------------------------------------------------------%
    % Panel (g)
    subplot(3,3,7);
    plot(y2t,'-k');
    title('(a) ARMA(2,0)');
    %ylabel('$y_t$');
    xlabel('$t$');
    axis tight
    box off;


    %--------------------------------------------------------%
    % Panel (h)
    subplot(3,3,8);
    bar(t,[1; acf2]);
    title('(b) ACF ARMA(2,0)');
    %ylabel('ACF');
    xlabel('lag');
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    axis tight
    box off;
     

    %--------------------------------------------------------%
    % Panel (i)
    subplot(3,3,9);    
    bar(t,[1; pacf2]);
    title('(b) PACF ARMA(2,0)');
    %ylabel('PACF');
    xlabel('lag');
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    axis tight
    box off;

    
    %laprint(1,'simrma','options','factory');

    %--------------------------------------------------------%
    %--------------------------------------------------------%
    
    figure(3)
    clf;
    
    % Panel (a)
    subplot(3,3,1);
    plot(y3t,'-k');
    title('(a) ARMA(0,1)');
    %ylabel('$y_t$');
    xlabel('$t$');
    axis tight
    box off;


    %--------------------------------------------------------%
    % Panel (b)
    subplot(3,3,2);
    bar(t,[1; acf3]);
    title('(b) ACF ARMA(0,1)');
    %ylabel('ACF');
    xlabel('lag');
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    axis tight
    box off;
     

    %--------------------------------------------------------%
    % Panel (c)
    subplot(3,3,3);    
    bar(t,[1; pacf3]);
    title('(b) PACF ARMA(0,1)');
    %ylabel('PACF');
    xlabel('lag');
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    axis tight
    box off;

  %--------------------------------------------------------%
  %--------------------------------------------------------%
    
    % Panel (d)
    subplot(3,3,4);
    plot(y4t,'-k');
    title('(a) ARMA(0,2)');
    %ylabel('$y_t$');
    xlabel('$t$');
    axis tight
    box off;


    %--------------------------------------------------------%
    % Panel (e)
    subplot(3,3,5);
    bar(t,[1; acf4]);
    title('(b) ACF ARMA(0,2)');
    %ylabel('ACF');
    xlabel('lag');
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    axis tight
    box off;
     

    %--------------------------------------------------------%
    % Panel (f)
    subplot(3,3,6);    
    bar(t,[1; pacf4]);
    title('(b) PACF ARMA(0,2)');
    %ylabel('PACF');
    xlabel('lag');
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    axis tight
    box off;

    %--------------------------------------------------------%
    %--------------------------------------------------------%
   
    % Panel (g)
    subplot(3,3,7);
    plot(y5t,'-k');
    title('(a) ARMA(1,1)');
    %ylabel('$y_t$');
    xlabel('$t$');
    box off;


    %--------------------------------------------------------%
    % Panel (h)
    subplot(3,3,8);
    bar(t,[1; acf5]);
    title('(b) ACF ARMA(1,1)');
    %ylabel('ACF');
    xlabel('lag');
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    axis tight
    box off;
     

    %--------------------------------------------------------%
    % Panel (i)
    subplot(3,3,9);    
    bar(t,[1; pacf5]);
    title('(b) PACF ARMA(1,1)');
    %ylabel('PACF');
    xlabel('lag');
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    axis tight
    box off;
   
   %laprint(1,'simarma','options','factory');


end
%
%------------------------- Functions -------------------------------------%
%
%-------------------------------------------------------------------------%
% ACF - based on running a sequence of OLS regressions 
%-------------------------------------------------------------------------%
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
%-------------------------------------------------------------------------%
% PACF - based on running a sequence of OLS regressions 
%-------------------------------------------------------------------------%
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

    