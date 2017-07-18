%=========================================================================
%
%   Program to simulate bilinear time series models
%
%=========================================================================

function nlm_bilinear( )

    clear all
    clc

    % Set seed of random number generator
    RandStream.setGlobalStream( RandStream('mt19937ar','seed',12) ); 

    % Parameters
    t = 200;

    % Simulate the models
    % Gamma = 0.0;
    y1 = bilinear(t,0.4,0.2,0.0);       
    
    % Gamma = 0.4
    y2 = bilinear(t,0.4,0.2,0.4);       
    
    % Gamma = 0.8
    y3 = bilinear(t,0.4,0.2,0.8);          

    % Gamma = 1.2
    y4 = bilinear(t,0.4,0.2,1.2); 
    
    
    %**********************************************************************
    %***
    %***     Generate graph
    %***
    %**********************************************************************

    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;
   
    s = 1:1:t;
    
    %--------------------------------------------------------%
    % Panel (a)
    subplot(2,2,1)
    plot(s,y1,'-k')
    title('(a) $\gamma = 0.0$');
    ylabel('$y$');
    xlabel('$t$');
    axis tight 
    box off

    %--------------------------------------------------------%
    % Panel (b)
    subplot(2,2,2)
    plot(s,y2,'-k')
    title('(b) $\gamma = 0.4$');
    ylabel('$y$');
    xlabel('$t$');
    axis tight 
    box off

%--------------------------------------------------------%
    % Panel (c)
    subplot(2,2,3)
    plot(s,y3,'-k')
    title('(c) $\gamma = 0.8$');
    ylabel('$y$');
    xlabel('$t$');
    axis tight 
    box off


%--------------------------------------------------------%
    % Panel (d)
    subplot(2,2,4)
    plot(s,y4,'-k')
    title('(d) $\gamma = 1.2$');
    ylabel('$y$');
    xlabel('$t$');
    axis tight 
    box off


end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Function to simulate a bilinear model
%-------------------------------------------------------------------------
function y = bilinear(t,phi,theta,gam)
  

    y = zeros(t+100,1) ;
    u = randn(t+100,1);
     
    for i = 2:t+100

        y(i) = phi*y(i-1) + u(i) + theta*u(i-1) + gam*y(i-1)*u(i-1);

    end
    y = y(101:end);

end

