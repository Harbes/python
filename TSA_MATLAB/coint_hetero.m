%=========================================================================
%
%   Program to investigate the effects of heteroskedasticity 
%   on the size of the trace test of cointegration
%
%=========================================================================
function coint_hetero( )

    clear all
    clc
      
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234) );

    reps  = 10000;
    n     = 4;
    model = 1;
    Tv    = [100,200,400,800,1600,3200];
    cv1   = [4.173,12.285,24.102,39.921,59.829,83.428];

    tr = zeros(reps,6); 

    for tc = 1:length(Tv) 
  
    t = Tv(tc);
  
    for rep = 1:reps 
        
        e = randn(t,n);

        % iid 
        tr(rep,1) = tracetest(cumsum(e),1,model);

        % GARCH(1,1)
        tr(rep,2) = tracetest(cumsum(garchdgp(e,0.3,0.6)),1,model);
        tr(rep,3) = tracetest(cumsum(garchdgp(e,0.3,0.69)),1,model);

        % Variance break 
        tr(rep,4) = tracetest(cumsum(varbreakdgp(e,0.1,2)),1,model);
        tr(rep,5) = tracetest(cumsum(varbreakdgp(e,0.5,2)),1,model);
        tr(rep,6) = tracetest(cumsum(varbreakdgp(e,0.9,2)),1,model);

    end
    format short g
    disp('          T           iid        GARCHa       GARCHb       Break-a      Break-b     Break-c');    
    disp( [t mean(tr>cv1(n)) ] );

    end

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Variance break DGP 
%-------------------------------------------------------------------------
function y = varbreakdgp(e,frac,newvar)

    [t,n] = size(e);
    tb    = floor(t*frac);
    sig2  = [ ones(tb,n); (newvar*ones(t-tb,n)) ];
    y     = e.*sig2;

end
%-------------------------------------------------------------------------
%  GARCH(1,1) DGP 
%-------------------------------------------------------------------------
function y = garchdgp(e,phi1,phi2)

    [t,n] = size(e);
    v2    =  e;
    sig2  =  ones(t,n);
    
    for j = 2:t;
        sig2(j,:) = (1-phi1-phi2) + phi1*v2(j-1,:).^2 + phi2*sig2(j-1,:);
        v2(j,:)   = e(j,:).*sqrt(sig2(j,:));
    end
    y = v2;
end
%-------------------------------------------------------------------------
%  Trace test
%-------------------------------------------------------------------------
function tr = tracetest(y,p,model)
    
        dy = trimr(y,1,0)-trimr(y,0,1); 
        z0 = trimr(dy,p-1,0);
        z1 = trimr(y,p-1,1);
        z2 = [];
        
        for j =1:p-1
        
            z2 = [ z2 trimr(dy,p-1-j,j)];
        end
    
        if model == 2
            
            z1 = [ trimr(y,p-1,1) ones(length(y)-p,1) ];
        
        elseif model == 3 
        
            z2 = [ones(length(y)-p,1) z2 ];   
    
        elseif model == 4
        
            z1 = [ z1 seqa(1,1,length(y)-p)'];
            z2 = [ones(length(y)-p,1) z2 ];   
        
        elseif model == 5; 
            
            z2 = [ones(length(y)-p,1) seqa(1,1,rows(y)-p)' z2];   
        end

    if size(z2,2) == 0 
        
        r0 = z0; 
        r1 = z1;
    
    else
        r0 = z0 - z2*(z2\z0); 
        r1 = z1 - z2*(z2\z1);
        
    end

    C     = inv(chol(r1'*r1));  
    lambda = eig(C'*r1'*r0*inv(r0'*r0)*r0'*r1*C);
    
    tr = -length(r0)*sum(log(1-flipud(lambda)));
end
