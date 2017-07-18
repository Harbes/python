% ========================================================================
%
%   Algorithm to estimate the FVAR model of Stock and Watson 
%   using yield data for the US.
%   US daily yields starting Jan 10th 1990 and ending May 31 2005
%
% ========================================================================

function lfac_fvar( )

    clear all
    clc
    
    % Load the data		
    
    % Yields 1yr, 2yr ... to 30yr    
    load('daily_finance.mat');
    
    yields = yields(:,1:30);
    
    % Standardize the data    
    y  = zscore(yields);
    [ aa,n ] = size( y );        % Number of dependent variables
    k = 3;                      % Number of factors (chosen to be 3 here) 

    lam = zeros( n,k );
    gam = zeros( n,n );
    
    lamnew = lam;
    gamnew = gam;

    maxiter = 100;       % Maximum number of iterations   
    tol     = 1e-3;      % Convergence criterion

    % Estimate F-VAR parameters using the Stock-Watson algorithm  

    iter = 1;
    
    while iter <= maxiter
        
        % Principal component analysis                                     
        [ va,ve ] = eigs(corrcoef( y ) );                              
        s         = y*va(:,1:k);
               
        % Regressions on measurement equations    
        for i = 1:n 
            ylag = trimr(y(:,i),0,1);
            b    =  [ trimr(s,1,0) ylag ]\trimr(y(:,i),1,0);
            lam(i,:) = b(1:k)';
            gam(i,i) = b(k+1)';
            y(:,i)   = y(:,i) - gam(i,i)*( [0.0 ; ylag]);              
        end
                
        % Check for convergence
        if sqrt(sum(vec(lam-lamnew).^2)) < tol && sqrt(sum(vec(gam-gamnew).^2)) < tol     

            break
        else
            lamnew = lam;
            gamnew = gam;
        end       
        iter = iter+1;  
    end

    % Regressions on factor equations 
    phi = trimr(s,0,1)\trimr(s,1,0);
    
    % Scale lam so shock is a positive shock to interest rates
    lam = -lam;
    
    disp( 'Convergence achieved after iteration');
    disp( iter );

    disp( 'Lambda');
    disp( lam );

    disp( 'Gamma' );
    disp( diag(gam) );

    disp( 'Phi');
    disp( phi );

    
    %**********************************************************************
    %***
    %***     Generate graphs
    %***
    %**********************************************************************

    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;

    t=1:30;
    ymax = 0.5;
    ymin = -0.5;
    
    %--------------------------------------------------------%
    
    plot(t,lam(:,1),'-.k',t,lam(:,2),'--k',t,lam(:,3),'-k','LineWidth',1);
    ylim([ymin ymax])
    %title('Three Factors')
    %legend('Level','Slope','Curvature','Location','SouthEast');
    legend boxoff;
    xlabel('Maturity');
    box off
    %axis tight
    % Print the tex file to the relevant directory
    %laprint(1,'figklevel','options','factory');


end