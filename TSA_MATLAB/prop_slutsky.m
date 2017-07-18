%=========================================================================
%
%   Program to demonstrate Slutsky's theorem by simulation 
%
%=========================================================================

clear all
clc

RandStream.setDefaultStream( RandStream('mt19937ar','seed',12) )

% Choose parameters
t  = 10;
n  = 50000;
mu = 2;

m1 = zeros(n,1);        
m2 = zeros(n,1);        

for i =  1:n

    % Generate exponential random numbers from a Gamma distribution
    y = mu*randg(1,[t 1]);         

    m1(i) = ( mean(y)/std(y) )^2;                                                     
    m2(i) = ( sqrt(t)*mean(y) )^2;                                                         
                                               

end
format ShortG

disp( [' Sample size    = ' num2str(t) ] );
disp( ' ' );
disp( ' Moment results for m1 ');
disp( ['Mean of m1     = ' num2str(mean(m1)) ] );
disp( ['Variance of m1 = ' num2str(std(m1)^2) ] );
disp( ' ' );
disp( ' Moment results for m2 ');
disp( ['Mean of m2     = ' num2str(mean(m2)) ] );
disp( ['Variance of m2 = ' num2str(std(m2)^2) ] );
