% ========================================================================
%
%      Program to generate GMM table and demonstrate empirical moments
%
% ========================================================================

clear all;
clc;


RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) );

t = 10;                         
%yt = round(5 + 2*randn(t,1));

yt = load('table.dat','-ascii'); 

m1 = mean(yt);                 

mt = [yt , (yt.^2) , (yt.^3) , (yt.^4) , (yt - m1) , (yt - m1).^2 , ...
      log(yt) , (1./yt)];

m = mean(mt);
disp([ mt;m ]);


disp([ 'Normal distribution (mu)    = ', num2str(m(1)) ]);
disp([ 'Normal distribution (sig2)  = ', num2str(m(6)) ]);

disp( ' ' );
disp([ 'Student t distribution (mu) = ', num2str(m(1)) ]);
disp([ 'Student t distribution (nu) = ', num2str(2*m(6)/(m(6)-1)) ]);


% Fourth central moment 
% (can be computed directly as c4 = mean((yt-m(1)).^4)
c4 = m(4)-4*m(3)*m(1)+6*m(2)*m(1)^2-4*m(1)*m(1)^3+m(1)^4;           

disp( ' ' );
disp([ 'Student t distribution (mue) = ', num2str(m(1)) ]);
% Solution based on solving a quadratic equation in nu
disp('Student t distribution (nu)' );                            
disp([ '     first root              = ', num2str((6+sqrt(36 - 4*(1-3/c4)*8))/(2*(1-3/c4))) ]);    
disp([ '     second root             = ', num2str((6-sqrt(36 - 4*(1-3/c4)*8))/(2*(1-3/c4))) ]);    

disp( ' ' );
disp([ 'Gamma distribution (alpha)   = ', num2str(m(1)/(m(1)-m(8))) ]);
disp([ 'Gamma distribution (beta)    = ', num2str((m(1)/(m(1)-m(8)))/m(1)) ]);

disp( ' ' );
disp([ 'Gamma distribution (alpha)   = ', num2str(m(1)^2/(m(2)-m(1)^2)) ]);
disp([ 'Gamma distribution (beta)    = ', num2str(m(1)/(m(2)-m(1)^2)) ]);

disp( ' ' );
disp([ 'Pareto distribution (alpha)  = ', num2str(m(1)/(m(1)-min(yt))) ]);

