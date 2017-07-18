%=========================================================================
%
%   The effect of a neglected break in trend on unit root tests 
%   DF and MZt tests with OLS and GLS de-trending
%
%=========================================================================
function unit_breakeffect( ) 

    clear all
    clc
    
    RandStream.setGlobalStream( RandStream('mt19937ar','seed',1234) )
    
    % Parameters
    t      = 200;
    reps   = 10000;
    tauv   = [0.25,0.5,0.75];
    breakv = [-0.4,-0.2,0.2,0.4];
    phiv   = [1,0.95,0.9,0.85];
    cv     = [-3.42 -2.86 -3.13 -2.86];
    
    rej  = zeros(length(phiv),4);
    s   = zeros(reps,4);
    x   = [ones(t,1) seqa(1,1,t)'];

    % Results with no break for comparison
    for phic = 1:length(phiv)
        for rep = 1:reps
            
            y    = recserar(randn(t,1),randn(1,1),phiv(phic));
            uols = glsdetrend(y,x,-t);
            ugls = glsdetrend(y,x,-13.5);
    
            s(rep,1) = trimr(adftests(uols,0),1,0);
            s(rep,2) = trimr(adftests(ugls,0),1,0);
            s(rep,3) = trimr(mtests(uols,0)',2,0);
            s(rep,4) = trimr(mtests(ugls,0)',2,0);
        end
    rej(phic,:) = mean(bsxfun(@lt,s,cv))';
    end
    disp('No Break')
    disp('    phi       DF-OLS    DF-GLS    MZt-OLS   MZt-GLS');
    disp([phiv' rej])
 

    % Result with break in dgp
    for tc = 1:length(tauv) 
        
        TB = floor(tauv(tc)*t);
        DT = [zeros(TB,1) ; seqa(1,1,t-TB)' ];
        x0 = [ x DT ];

        for bc = 1:length(breakv) 
     
            beta = [ 0; 0; breakv(bc) ];
            for phic = 1:length(phiv) 
      
                for rep = 1:reps 
                    
                    u = recserar(randn(t,1),randn(1,1),phiv(phic));
                    y = x0*beta+u;
                    uols = glsdetrend(y,x,-t);
                    ugls = glsdetrend(y,x,-13.5);
            
                    s(rep,1) = trimr(adftests(uols,0),1,0);
                    s(rep,2) = trimr(adftests(ugls,0),1,0);
                    s(rep,3) = trimr(mtests(uols,0)',2,0);
                    s(rep,4) = trimr(mtests(ugls,0)',2,0);
                end
                rej(phic,:) = mean(bsxfun(@lt,s,cv))';
            end
        end
    disp( ' ' )
    disp(['Break point    = ', num2str(tauv(tc)) ])
    disp(['Break size     = ', num2str(breakv(bc)) ])
    disp('    phi       DF-OLS    DF-GLS    MZt-OLS   MZt-GLS');
    disp([phiv' rej])

    end
end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Detrending function: 
%       cbar = -7 constant; 
%       cbar = -13.5 linear trend  
%       cbar = -T for OLS detrending
%-------------------------------------------------------------------------
function [ u,b ] = glsdetrend( y,x,cbar )

    t = length( y );
    yc = [ y(1)  ; (trimr(y,1,0)-(1+cbar/t)*trimr(y,0,1)) ];
    xc = [x(1,:) ; (trimr(x,1,0)-(1+cbar/t)*trimr(x,0,1)) ];
    
    b  = xc\yc;
    u  = y - x*b;
 
end
%-------------------------------------------------------------------------
%  ADF coefficient and t tests: u must be already de-trended
%-------------------------------------------------------------------------
function adfstats= adftests(u,k)


    du = trimr(u,1,0) - trimr(u,0,1);
    x  = trimr(u,0,1);

    % Set up lags of du
    if k>0;
    
        ldu = lagmatrix(du,1:1:k);
        x   = trimr( [x ldu],k,0 );
    end
    
    xxi = inv(x'*x); 
    b   = xxi*x'*trimr(du,k,0); 
    e   = trimr(du,k,0)-x*b; 
    s2  = e'*e/length(e);

    adfstats = [length(u)*b(1) ; b(1)/sqrt(s2*xxi(1,1))];

end

%-------------------------------------------------------------------------
%  M tests
%-------------------------------------------------------------------------
function tests = mtests(u,k)

    s2  = ar1rvar(u,k);
    n   = length(u);
    tmp = sum(u(1:n-1).^2);
    
    %disp(['Sum of u(1) to u(t-1) = ',num2str(tmp) ])
    %disp(['Last value: u(t)      = ',num2str(u(n)) ])
    %disp(' ') 

    u2  = tmp/n^2;     
    mza = (u(n)^2/n-s2)/(2*u2);     
    msb = sqrt(u2/s2);
    mzt = msb*mza;
    
    tests = [mza msb mzt];
 
end
%-------------------------------------------------------------------------
%  Autoregressive long run variance estimator
%-------------------------------------------------------------------------
function s2 = ar1rvar(u,k)

    du = trimr(u,1,0)-trimr(u,0,1); 
    x  = trimr(u,0,1);

     if k > 0 
        x  = [ x lagmatrix(du,seqa(1,1,k)) ];
        x(any(isnan(x),2),:) = [];
     end


    b = x\trimr(du,k,0); 
    e = trimr(du,k,0) - x*b; 
    s2 = e'*e/length(e);

    if k>0 
        s2 = s2/(1-sum(trimr(b,1,0)))^2; 
    end

end

