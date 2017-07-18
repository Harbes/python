%=========================================================================
%
%   KPSS asymptotic critical values
%
%=========================================================================
function unit_kpss_cv(  )

    clear all
    clc

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234) )
    
    reps = 100000;
    t    = 1000;
    
    % Constant
    cv = kpsscv(ones(t,1),t,reps);
    disp(' Constant ' )
    disp('    0.10      0.05      0.01 ')
    disp('    -------------------------')
    disp( cv)
    
    % Constant and trend
    cv = kpsscv([ones(t,1) seqa(1,1,t)'],t,reps);   
    disp( ' ' )
    disp(' Constant and trend' )
    disp('    0.10      0.05      0.01 ')
    disp('    -------------------------')
    disp( cv)

      
    tau=0.25; 
    cv = kpsscv([ones(t,1) (seqa(1,1,t)>floor(tau*t))'],t,reps);
    disp(' tau = 0.25 : Constant and trend' )
    disp('    0.10      0.05      0.01 ')
    disp('    -------------------------')
    disp( cv)

    tau=0.5; 
    cv = kpsscv([ones(t,1) (seqa(1,1,t)>floor(tau*t))'],t,reps);
    disp(' tau = 0.5 : Constant and trend' )
    disp('    0.10      0.05      0.01 ')
    disp('    -------------------------')
    disp( cv)
   
    tau=0.75; 
    cv = kpsscv([ones(t,1) (seqa(1,1,t)>floor(tau*t))'],t,reps);
    disp(' tau = 0.75 : Constant and trend' )
    disp('    0.10      0.05      0.01 ')
    disp('    -------------------------')
    disp( cv)

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  KPSS test
%-------------------------------------------------------------------------
function cv = kpsscv(x,t,reps)

    kpss = zeros(reps,1);
    for j = 1:reps
        
        y = randn(t,1);
        z = y-x*(x\y);
        S = cumsum(z);
        
        kpss(j) = (S'*S)/(z'*z*t);
    end
    
    cv = quantile(kpss,[0.90 0.95 0.99]);

end

