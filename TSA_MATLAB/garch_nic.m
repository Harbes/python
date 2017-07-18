%=========================================================================
%
%      Program to compute and plot the News Impact Curve 
%
%=========================================================================

function garch_nic( )

    clear all
    clc
        
    % Set random number generator
    RandStream.setGlobalStream( RandStream('mt19937ar','seed',1) );
  
    % Parameter values
    a0 = 0.1;
    a1 = 0.5;
    b1 = 0.4;


    % Compute news impact curve of an ARCH(1) model
    u = -5:0.1:5;

    a0 = 1;
    a1 = 0.2;
    sig21 = a0 + a1*u.^2;

    a0 = 1;
    a1 = 0.5;
    sig22 = a0 + a1*u.^2;

    a0 = 1;
    a1 = 0.8;
    sig23 = a0 + a1*u.^2;

   
    %*********************************************************************
    %**     Generate graph of news impact curve
    %*********************************************************************
    
    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    fig1=figure(1);
    clf;


    plot( u,sig21,'-k',u,sig22,'--k',u,sig23,':k','LineWidth',0.75);
    ylabel('$h_t$');
    xlabel('$u$');
    legend( '$\alpha=0.2$','$\alpha=0.5$','$\alpha=0.8$' );  
    legend boxoff
    legend2latex(fig1); 
    %laprint(1,'nic','options','factory');



end