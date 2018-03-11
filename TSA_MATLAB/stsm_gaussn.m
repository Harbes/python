%=========================================================================
%
%   Program to estimate an ARMA(1,1) using the GAUSS-NEWTON algorithm
%
%=========================================================================
clear all;
clc;

%RandStream.setDefaultStream( RandStream('mt19937ar','seed',12345) );

% Alternative starting values
theta = [0.2 ; 0.2 ; 0.2 ];       
%theta = [0.0 ; 0.0 ; 0.0 ];       
%theta = [0.0 ; 0.3 ; 1.1 ]; 
%theta = [0.0 ; 1.1 ; 0.3 ];

% Alternative sample sizes
t = 200;
%t = 500;

% Simulate the data (discard first 100 simulated observations)
mu   = 0.0;            
phi1 = 0.8;           
si1  = 0.3;
vt   = randn(t+101,1);
yt   = recserar( mu + trimr(vt,1,0) + si1*trimr(vt,0,1) , 0.0 , phi1);  
yt   = trimr(yt,100,0);   

% Gauss Newton algorithm
crit  = 0.00001;                % Convergence criterion          
maxit = 20;                      % Maximum number of iterations    


i = 1;
while i <= maxit
    
    % vt is being redefined in the gauss-newton algorithm   
    vt  = recserar( trimr(yt,1,0) - theta(1) - theta(2)*trimr(yt,0,1) , 0.0 , -theta(3));        
    z1t = recserar( ones(length(vt),1) , 0.0 , -theta(3) );
    z2t = recserar( trimr(yt,0,1) , 0.0 , -theta(3) );
    z3t = recserar( [ 0.0 ; trimr(vt,0,1) ] , 0.0 , -theta(3));
    zt  = [ z1t  z2t  z3t ];

    % Remove all starting zeros   
    vt  = trimr(vt,2,0);       
    zt  = trimr(zt,2,0);

 
    disp(['Iteration        = ', num2str(i) ]);
    disp(['Parameter vector = ', num2str(theta') ]);
    disp(' ');

    dtheta = zt\vt;

    % Check convergence and update parameters
    if dtheta'*dtheta < crit
        break
    else
        theta   = theta  + dtheta;  
        if i == maxit
            disp(['Failed to converge after iteration   ', num2str(maxit) ]);
        end

    end
    i = i + 1;
end

sig2 = (vt'*vt)/t;
vc   = sig2*inv(zt'*zt);
se   = sqrt(diag(vc));
disp(' ')
disp('Parameters Std. errors   t-stats')
disp( [theta  se  theta./se ]);
disp(' ')
disp('Estimated asymptotic covariance matrix')
disp( vc )

% Wald test    
w = theta(3)^2/vc(3,3);

disp(' ')
disp(['Wald statistic (si1 = 0.0)     = ', num2str(w) ]);
disp(['p-value                        = ', num2str(1-chi2cdf(w,1)) ]);
disp(' ')

% LM test   
% First regression
y = trimr(yt,1,0);   
x = [ones(length(y),1)   trimr(yt,0,1) ];
e = y - x*(x\y);

% Second regression
y   = trimr(e,1,0);    
x   = [ones(length(y),1)   trimr(yt,1,1)   trimr(e,0,1) ];
e   = y - x*(x\y);
y   = bsxfun(@minus,y,mean(y));
rsq = 1 - (e'*e)/(y'*y);     
lm  = length(yt)*rsq;

disp(' ')
disp(['LM statistic (si1 = 0.0)     = ', num2str(lm) ]);
disp(['p-value                        = ', num2str(1-chi2cdf(lm,1)) ]);
disp(' ')

