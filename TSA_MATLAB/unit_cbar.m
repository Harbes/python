%=========================================================================
%
%   Program to find the value of cbar such that the power of the point
%   optimal unit root test against phi = 1+cbar/T has asymptotic power of 0.5
%
%=========================================================================
function unit_cbar( )

    clear all
    clc
    
    reps = 1000;
    t    = 100;

    cbarc = cbar(ones(t,1),-6,reps,t);
    cbart = cbar([ones(t,1) seqa(1,1,t)'],-12,reps,t);

    tauv    = seqa(0.15,0.05,15);
    cbartau = zeros(length(tauv),1);

    for tc = 1:length(tauv)
    
        TB = floor(tauv(tc)*t);
        cbartau(tc) = cbar([ones(t,1) seqa(1,1,t)' [zeros(TB,1);seqa(1,1,t-TB)']],-13,reps,t);
    end
            
    disp(['Constant               = ',num2str(cbarc) ])
   	disp(['Linear trend           = ',num2str(cbart) ])
    disp( ' ' );
    disp('Break in trend')
    disp('tau        cbar')
    disp([tauv' cbartau])

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Search for cbar for a given x
%-------------------------------------------------------------------------
function cbarc = cbar(x,cbar0,reps,t)

    cbarv = cbar0; 
    powerv = power(x,cbarv,reps,t);
    
    while powerv(1) < 0.5

        powerv = [ power(x,cbarv(1),reps,t) ; powerv ];
        cbarv  = [ cbarv(1)-0.5  ;  cbarv];
        
    end
    
    cbarc = cbarv(2) - 0.5*(0.5-powerv(2)/(powerv(1)-powerv(2)));

end
%-------------------------------------------------------------------------
%  Power envelope
%-------------------------------------------------------------------------
function  pow = power(x,c,reps,t)

    RandStream.setGlobalStream( RandStream('mt19937ar','seed',42) )       

    p0 = zeros(reps,1); 
    p1 = zeros(reps,1);

    for rep = 1:reps 
        u  = randn(t,1);
        y0 = cumsum(u);
        y1 = recserar(u,u(1),1+c/t);
        p0(rep) = Ptest(y0,x,c,0);
        p1(rep) = Ptest(y1,x,c,0);
    end
    
    pow = mean(p1 < quantile(p0,0.05));
end
%-------------------------------------------------------------------------
%  P-test
%-------------------------------------------------------------------------
function pt = Ptest(y,x,cbar,k)

    n  = length(y);
    uc = glsdetrend(y,x,cbar); 
    u0 = glsdetrend(y,x,0);
    s2 = ar1rvar(uc,k);
    uc = [ uc(1) ; trimr(uc,1,0)-(1+cbar/n)*trimr(uc,0,1) ];
    u0 = [ u0(1) ; trimr(u0,1,0)-trimr(u0,0,1) ];

    pt = (uc'*uc-(1+cbar/n)*u0'*u0)/s2;

end
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
%  ar1rvar
%-------------------------------------------------------------------------
function s2 = ar1rvar(u,k)

    du = trimr(u,1,0)-trimr(u,0,1); 
    x  = trimr(u,0,1);

    if k>0 
        
        tmp = [x lagmatrix(du,1:1:k)];
        
        x = trimr( tmp,k,0 ); 
    end

    b = x\trimr(du,k,0); 
    e = trimr(du,k,0) - x*b; 
    s2 = e'*e/length(e);

    if k>0
        
        s2 = s2/(1-sum(trimr(b,1,0)))^2; 
    
    end
    

end
