%=========================================================================
%
%   Program to simulate the OLS estimator of the AR(1) model
%
%=========================================================================

function nts_distribution ( )

    clear all;
    clc;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',5) );

    nreps = 50000;           
    t     = [ 50,100,200,400,800,1600 ];
   
    % Compute statistics on the sampling distribution (stationary)
    disp(['Statistics of the sampling distribution for phi = ',num2str(0.8)]);
    
    xi = (-5:0.05:5)';
    fs  = zeros( length(xi),length(t) );
       
    for j = 1:length(t)
              
        zs   = ar1sim(0.8,t(j),nreps);  
        bias = mean( zs  ); 
        stdv = std( zs  ); 
        skew = mean( zs.^3 )/stdv^3; 
        kurt = mean( zs.^4 )/stdv^4;
     
        disp(['t                         = ',num2str( t(j) )]);
        disp(['Bias                      = ',num2str( bias )]);
        disp(['Standard deviation        = ',num2str( stdv )]);
        disp(['Skewness                  = ',num2str( skew )]);
        disp(['Kurtosis                  = ',num2str( kurt )]);
        %disp(['Proportion of stats < 0.0 = ',num2str( mean(z<=0) )]);
        disp(' ' );
            
        fs(:,j) = ksdensity( zs,xi );
            
    end
        
    tmps = normpdf( xi );

    % Compute statistics on the sampling distribution (stationary)
    disp(['Statistics of the sampling distribution for phi = ',num2str(1.0)]);
    
    xin = (-15:0.01:5)';
    fn  = zeros( length(xin),length(t) );
    
    for j = 1:length(t)
        
        zn   = ar1sim(1.0,t(j),nreps);      
        bias = mean( zn ); 
        stdv = std( zn ); 
        skew = mean( zn.^3 )/stdv^3; 
        kurt = mean( zn.^4 )/stdv^4;
     
        disp(['T                         = ',num2str( t(j) )]);
        disp(['Bias                      = ',num2str( bias )]);
        disp(['Standard deviation        = ',num2str( stdv )]);
        disp(['Skewness                  = ',num2str( skew )]);
        disp(['Kurtosis                  = ',num2str( kurt )]);
        %disp(['Proportion of stats < 0.0 = ',num2str( mean(z<=0) )]);
        disp(' ' );
            
        fn(:,j) = ksdensity( zn,xin );
        
    end    
    tmpn = normpdf( xin );
    
    %**********************************************************************
    %***
    %***     Generate graph
    %***
    %**********************************************************************

    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;
    
    %--------------------------------------------------------%
    % Panel (a)
    subplot(1,2,1)
    plot(xi,tmps,'-k',    ...
         xi,fs(:,1),':k', ...
         xi,fs(:,3),'--k');
    title('(a) $\phi=0.8$');
    xlabel('$z$')
    ylabel('$f(z)$')
    set(gca,'YTick',[0.1 0.2 0.3 0.4])
    h1 = subplot(121);
    axis(h1,'tight')
    axis
    box off
      
    %--------------------------------------------------------%
    % Panel (b)
    subplot(1,2,2)
    plot(xin,tmpn,'-k',    ...
         xin,fn(:,1),':k', ...
         xin,fn(:,3),'--k');
    title('(b) $\phi=1.0$');
    xlabel('$z$')
    ylabel('$f(z)$')
    set(gca,'YTick',[0.1 0.2 0.3 0.4])
    box off
  
    % Print the tex file to the relevant directory
    %laprint(1,'AR1dist','options','factory');
    
end

%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
% AR1 Procedure  
%-------------------------------------------------------------------------
function z = ar1sim( phi,t,nreps)
   
    phihat = zeros(nreps,1);

    %  Loop over reps simulating and estimating AR(1) model
    for j = 1:nreps
        
        v         = randn(t+1,1);
        y         = trimr(recserar(v,0,phi),1,0);   
        phihat(j) = y(1:t-1)\y(2:t);  
        
    end

    if phi<1;
        z = sqrt(t)*(phihat - phi)/sqrt(1 - phi^2); % Stationary case
    else
        z = t*(phihat - phi);                       % Nonstationary case
    end

end   
 
