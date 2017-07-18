%=========================================================================
%
%     Plot a near unit root process
%
%=========================================================================
function unit_nearunitroot( )

    clear all;
    clc;
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) );

    t    = 200;                            
    sig2 = 0.1;
    
    % Simulate AR1 with phi=1.0 
    mue  = 0.0;  
    phi1 = 1.0;  
    phi2 = 0.0; 
    si1  = 0.0; 
    si2  = 0.0;

    vt  = sqrt(sig2)*randn(t+101,1);      
    yt  = recserar( mue + trimr(vt,2,0) + si1*trimr(vt,1,1) + si2*trimr(vt,0,2) , [0.0;0.0] , [phi1; phi2] );   
    y1t = trimr(yt,100,0);                               

    % Simulate AR1 with phi=0.99      
    mue  = 0.0;  
    phi1 = 0.99;  
    phi2 = 0.0; 
    si1 = 0.0; 
    si2 = 0.0;

    vt  = sqrt(sig2)*randn(t+101,1);     
    yt  = recserar( mue + trimr(vt,2,0) + si1*trimr(vt,1,1) + si2*trimr(vt,0,2) , [0.0;0.0] , [phi1;phi2] );   
    y2t = trimr(yt,100,0);                               


    %**********************************************************************
    %***
    %***     Generate graph
    %***
    %**********************************************************************

    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;
    
    x = 1:1:length(y1t);
    %--------------------------------------------------------%
    % Panel (a)
    subplot(1,2,1)
    plot( x,y1t,'-k' );
    title('(a) $\phi=1.0$')
    xlabel('$t$');
    box off
    ylim([-3 3]);
    

    x = 1:1:length(y2t);
    %--------------------------------------------------------%
    % Panel (b)
    subplot(1,2,2);
    plot( x,y2t,'-k' );
    title('(b) $\phi=0.99$')
    xlabel('$t$');
    box off
    ylim([-3 3]);
    


    %laprint(1,'nearunitroot','options','factory'); 
    
end

