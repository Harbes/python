%========================================================================
%
%      Program to estimate a exponential model 
%
%========================================================================

clear all
clc

% Data
y = [ 3.5, 1.0, 1.5 ]';
t = length(y);       	

niter = 10;  
% Newton-Raphson algorithm  
disp('Newton Raphson algorithm');
disp(' ');

theta = 1.0;
for k = 1:niter
    
    g     = -1/theta + mean(y)/theta^2;
    h     = 1/theta^2 - 2*mean(y)/theta^3;
    lnl   = -log(theta) - mean(y)/theta;
    theta = theta - inv(h)*g;
 
    disp( '     Iter      G         J        logL      th(k)') ;
    disp( [k g  g'*g lnl theta ] );

end

disp(' ');

% Scoring algorithm 
disp('Scoring algorithm');
disp(' ');

theta = 1.0;
for k = 1:niter
    
    g     = -1/theta + mean(y)/theta^2;
    i     = 1/theta^2;
    lnl   = -log(theta) - mean(y)/theta;
    theta = theta + inv(i)*g;
    
    disp( '     Iter      G         J        logL      th(k)') ;
    disp( [k g g'*g lnl theta ] );

end

disp(' ');

% BHHH algorithm 
disp('BHHH algorithm');
disp(' ');

theta = 1.0;
for k = 1:niter
    
    gt    = -1/theta + y/theta^2;
    g     = mean(gt);
    j     = gt'*gt/t;
    lnl   = -log(theta) - mean(y)/theta;    
    theta = theta + inv(j)*g;
    
    disp( '     Iter      G         J        logL      th(k)') ;
    disp( [k g j lnl  theta ] );

end

% BHHH algorithm with squeezing  
disp('BHHH algorithm with squeezing');
disp(' ');

theta = 1.0;               

for k = 1:niter

    gt  = -1/theta + y/theta^2;
    g   = mean(gt);
    j   = gt'*gt/t; 
    lnl = -log(theta) - mean(y)/theta;

    thetaold = theta;
    lnlold   = lnl;
    theta    = theta + inv(j)*g;
    lnl      = -log(theta) - mean(y)/theta;
    
    if lnl > lnlold;      
        disp( ['Full iteration step successful at iteration k =' num2str(k) ]);
    else
        disp( ['Full iteration step not successful at iteration k =' num2str(k) ]);

                for m = 2:10

                    lambda = 1/m;
                    disp( ['Squeezing with step = ' num2str(lambda) ]);
                    theta = thetaold + lambda*inv(j)*g;
                    lnl   = -log(theta) - mean(y)/theta;
    
                    if lnl > lnlold
                        break                      
                    end                  
                end
    end
    disp( '     Iter     G         J        logL      th(k)') ;
    disp( [k g j lnl theta ] );


end
% BFGS algorithm   
disp('BFGS algorithm');
disp(' ');

theta = 1.5;               
h     = -1;

niter = 9;
for k = 1:niter

    g     = -1/theta + mean(y)/theta^2;     
    lnl   = -log(theta) - mean(y)/theta;
    theta = theta - inv(h)*g;

    disp( '     Iter       g         h      logL      th(k)') ;
    disp( [k g h lnl theta ] );
 
    dtheta = -inv(h)*g;    
    gold   = g;
    g      = -1/theta + mean(y)/theta^2;
    dg     = g - gold;         
    h      = h - ( h*dtheta*dg' + dg*dtheta'*h )/(dg'*dtheta) ...
             + ( 1 + (dtheta'*h*dtheta)/(dg'*dtheta) )*(dg*dg')/(dg'*dtheta);

end

 