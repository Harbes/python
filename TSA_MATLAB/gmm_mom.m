%=========================================================================
%
%   Program to compute method of moment estimates of various distributions
%
%=========================================================================
clear all
clc

RandStream.setDefaultStream( RandStream('mt19937ar','seed',12) );

t = 10;                            

% Generate data and compute sam
yt = round(5 + 2*randn(t,1));    
m1 = mean(yt);  
mt = [yt , (yt.^2) , (yt.^3) , (yt.^4) , (yt - m1) , (yt - m1).^2 , log(yt) , (1./yt)];
m  = mean(mt); 

disp( [mt; m]);
disp( ' ' );

% Compute various moments used for method of moments estimation
m2 = mean(yt.^2);
c2 = mean((yt - m1).^2);
c4 = m(4)-4*m(3)*m(1)+6*m(2)*m(1)^2-4*m(1)*m(1)^3+m(1)^4;      
h1 = mean(1./yt);

disp(['Normal distribution (mu)     = ',num2str(m1) ]);
disp(['Normal distribution (sig2)   = ',num2str(c2) ]);
disp( ' ' );

disp(['Normal distribution (mu)    = ',num2str(m1) ]);
disp(['Normal distribution (sig2)  = ',num2str(sqrt(c4/3)) ]);
disp( ' ' );

disp(['Student t distribution (mu) =  ',num2str(m1) ]);
disp(['Student t distribution (nu) = ',num2str(2*c2/(c2-1)) ]);
disp( ' ' );

disp(['Student t distribution (mue) = ',num2str(m1) ]);
disp('Student t distribution (nu) based on solving a quadratic equation in nu')
disp(['     first root              = ',num2str((6+sqrt(36 - 4*(1-3/c4)*8))/(2*(1-3/c4))) ]);    
disp(['     second root             = ',num2str((6-sqrt(36 - 4*(1-3/c4)*8))/(2*(1-3/c4))) ]);    
disp( ' ' );

disp(['Gamma distribution (alpha)   = ',num2str(m1/(m1-h1)) ]);
disp(['Gamma distribution (beta)    = ',num2str((m1/(m1-h1))/m1) ]);
disp( ' ' );

disp(['Gamma distribution (alpha)   = ',num2str(m1^2/(m2-m1^2)) ]);
disp(['Gamma distribution (beta)    = ',num2str(m1/(m2-m1^2)) ]);
disp( ' ' );

ymin = min(yt);
disp(['Pareto distribution (alpha)  = ',num2str(m1/(m1-ymin)) ]);
disp( ' ' );

% Solve a cubic equation to find the roots of alpha
alpha_roots = roots( [1 ; -4 ; 5-ymin^2/m2 ; -2] );        

 % The first root is real whereas the other two roots are the complex conjugates
disp(['Pareto distribution (alpha)  = ',num2str(alpha_roots(1)) ]);     



