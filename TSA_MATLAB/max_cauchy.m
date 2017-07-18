%=========================================================================
%
%   Program to estimate a Cauchy model (one iteration)
%   using the NEWTON-RAPHSON, SCORING and BHHH algorithms
%
%=========================================================================

clear all
clc

% Load data                     
yt = [ 2, 5, -2, 3, 3 ];
t  = length(yt);      


theta0 = 3;        


% Newton-Raphson
g = 2*sum( (yt - theta0)./(1 + (yt - theta0).^2) );
h = 2*sum( ((yt - theta0).^2 - 1)./(1 + (yt - theta0).^2).^2 );

thetaNR = theta0 - inv(h)*g;


% Method of Scoring
g = 2*sum( (yt - theta0)./(1 + (yt - theta0).^2) );
i = t/2;

thetaSC = theta0 + inv(i)*g;


% BHHH
gt = 2*( (yt - theta0)./(1 + (yt - theta0).^2) );

thetaBH = theta0 + inv(gt*gt')*sum(gt);


disp( ['Newton-Raphson     ' num2str(thetaNR) ] );
disp( ['Method of Scoring  ' num2str(thetaSC) ] );
disp( ['BHHH               ' num2str(thetaBH) ] );


disp( ' ' );
disp( ['Gradient                       ' num2str(g) ]);
disp( ['Standard error (Hessian)       ' num2str(sqrt(-inv(h))) ] );
disp( ['Standard error (Information)   ' num2str(sqrt(inv(i))) ] );
disp( ['Standard error (OPG)           ' num2str(sqrt(inv(gt*gt'))) ] );
disp( ' ' );



