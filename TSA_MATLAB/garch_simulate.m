%=========================================================================
%
%      Program to simulate a garch model 
%
%=========================================================================

function garch_simulate( )

    clear all
    clc
        
    % Set random number generator
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',12345) );
     
    nobs = 100000;
    mlag = 20;

    z  = randn( nobs+1000,1 );        
    u  = z;                         
    y  = z;
    h  = ones( nobs+1000,1 );

    % Choose parameter values  
    %a0 = 0.1; a1 = 0.7; b1 = 0.2;
    %a0 = 0.05; a1 = 0.15; b1 = 0.8;
    a0 = 0.05; a1 = 0.05; b1 = 0.9;

    h(1) = a0/(1-a1-b1);
    u(1) = normrnd(0,a0/(1-a1-b1));
    % Generate data 
    for t = 2:nobs+1000

        h(t) = a0 + a1*u(t-1)^2 + b1*h(t-1);    % Conditional variance  
        u(t) = z(t)*sqrt( h(t) );               % Disturbance term     
        y(t) = u(t);                            % Returns
    end
    y = trimr( y,1000,0 );  
    y = y - mean(y);

    
    %load test;
    % Compute the autocorrelation of returns squared
    acfy1    = zeros( mlag+1,1 );
 
     
    for j = 0:mlag 
    
        acfy1(j+1)  = acfunction( y.^2,j );
    end
 

    % Parameter values second model
    a0 = 0.05;
    a1 = 0.15;
    b1 = 0.80;
    
    
    % Generate data 
    for t = 2:nobs+1000

        h(t) = a0 + a1*u(t-1)^2 + b1*h(t-1);    % Conditional variance  
        u(t) = z(t)*sqrt( h(t) );               % Disturbance term     
        y(t) = u(t);                            % Returns
    end
    y = trimr( y,1000,0 );  
    y = y - mean(y);

    % Compute the autocorrelation of returns squared
    acfy2    = zeros( mlag+1,1 );    
    for j = 0:mlag 
    
        acfy2(j+1)  = acfunction( y.^2,j );
    end
    
    %*********************************************************************
    %**     Generate graphs
    %*********************************************************************
    
    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;
  
    lags = 0:1:mlag;
    
    subplot(2,1,1)
    plot(lags,[ acfy1 zeros(mlag+1,1) ] ,'-k','LineWidth',0.75);
    ylim([ -0.5 1 ] );
    xlim( [0 10] );
    title('(a)');
    ylabel('ACF $y_{1t}$');
    xlabel('$p$');
    box off;
    %axis tight;

    subplot(2,1,2)
    plot(lags,[ acfy2 zeros(mlag+1,1) ],'-k','LineWidth',0.75);
    ylim( [-0.5 1 ] );
    xlim( [0 10] );
    title('(b)');
    ylabel('ACF $y_{2t}$');
    xlabel('$p$');
    box off;
    %axis tight;

    %laprint(1,'volsim','options','factory');

end

%-------------------------------------------------------------------------
%   Subroutine: Compute the ACF at kag p
%-------------------------------------------------------------------------

function r = acfunction( x,p )

    d = x - mean( x );
    r = sum( trimr(d,p,0).*trimr(d,0,p) ) / sum(d.^2);

end






