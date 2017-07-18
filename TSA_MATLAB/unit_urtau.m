%=========================================================================
%
%   Compute tau for Union of Rejections MZt unit root test
%
%=========================================================================
function unit_urtau( )

    clear all
    clc

	RandStream.setGlobalStream( RandStream('mt19937ar','seed',1234) );
 
    t    = 1000; 
    reps = 100000;

    % Constant case
    x     = ones(t,1); 
    cbar  = -7; 
    cvols = -2.467; 
    cvgls = -1.944;

    % Constant and linear trend
    % x     = [ones(t,1) seqa(1,1,t)];
    % cbar  = -13.5; 
    % cvols = -3.130; 
    % cvgls = -2.867;


    % Obtain replications of MZt-OLS and MZt-GLS tests
    MZols = zeros(reps,1); 
    MZgls = zeros(reps,1);

    for rep = 1:reps
  
        disp('Obtaining replications of MZt-OLS and MZt-GLS tests ...')
        y          = cumsum(randn(t,1));
        MZols(rep) = trimr(mtests(glsdetrend(y,x,-t),0)',2,0);
        MZgls(rep) = trimr(mtests(glsdetrend(y,x,cbar),0)',2,0);
    end

    % Search for tau such that UR test has size of 0.05
    rejtau = 1; 
    tau    = 1;
    
    % Second decimal place:
    while rejtau >= 0.05
    
        rejtau = mean(MZols < (tau*cvols) | MZgls < (tau*cvgls));
        tau    = tau + 0.01;
    end

    rejtau = 1; 
    tau = tau - 0.02;
    
    % Third decimal place:
    while rejtau >= 0.05
        
        rejtau = mean(MZols < (tau*cvols) | MZgls < (tau*cvgls));
        tau    = tau + 0.001;
    end

    disp(['tau     = ',num2str(tau-0.001) ])

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
%  M tests
%-------------------------------------------------------------------------
function tests = mtests(u,k)

    s2  = ar1rvar(u,k);
    n   = length(u);
    tmp = sum(u(1:n-1).^2);
    
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