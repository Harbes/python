%=========================================================================
%
%   Program to demonstrate the consistency of GMM
%
%=========================================================================

function gmm_consistency( )

    clear all;
    clc;
    format short;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234) )


    % Population distribution: T=100000  
    t      = 100000;
    alpha0 = 10.0;
    alpha  = 6:0.1:14;
    y      = gamrnd(alpha0,1,t,1);  
    Q0     = zeros( length(alpha),1 );
    
    for k = 1:length(alpha)
          
        Q0(k) = gmmcrit( alpha(k),y );
        
    end;

    
    
    % Finite sample distributions
    Q = zeros( length(alpha),4 );
    t = 50;

     for k = 1:length(alpha)
        
        tmp    = y(1:t);  
        Q(k,1) = gmmcrit( alpha(k),tmp );
        
     end

    t = 100;
    for k = 1:length(alpha)
        
        tmp  = y(1:t);  
        Q(k,2) = gmmcrit( alpha(k),tmp );
        
    end
    
    t = 200;
    for k = 1:length(alpha)
        
        tmp  = y(1:t);  
        Q(k,3) = gmmcrit( alpha(k),tmp );
        
    end
    
    t = 400;
    for k = 1:length(alpha)
        
        tmp  = y(1:t);  
        Q(k,4) = gmmcrit( alpha(k),tmp );
        
    end

    %**********************************************************************
    %***
    %***     Generate graphs
    %***
    %**********************************************************************

    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;

    plot(alpha,Q0,'-k', ...
         alpha,Q(:,1),'-.k', ...
         alpha,Q(:,2),'--k');
    ylabel('$Q(\alpha)$');
    xlabel('$\alpha$');
    box off;
    
    laprint(1,'gmmconsistency','options','factory');

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
% GMM objective function  
%-------------------------------------------------------------------------
function Q = gmmcrit( alpha,y )

    d1 = y - alpha;
    d2 = y.^2 - alpha*(alpha+1);
    d3 = 1./y - 1/(alpha-1);
     
    d = [d1 d2 d3];
    g = mean(d);
    w = d'*d/length(d);

    Q = g*inv(w)*g';

end



