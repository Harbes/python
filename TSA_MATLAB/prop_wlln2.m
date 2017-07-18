%=========================================================================
%
%   Program to demonstrate by simulation the necessary and sufficient 
%   conditions of the weak law of large numbers
%
%=========================================================================

clear all
clc

RandStream.setDefaultStream( RandStream('mt19937ar','seed',12) )

% Choose parameters
t = 100;
n = 50000;

m1 = zeros(n,1);        
m2 = zeros(n,1);        
m3 = zeros(n,1);       
m4 = zeros(n,1);       


for i =  1:n

    y = rand(t,1) - 0.5;        % Simulate y from a uniform distribution (-0.5,0.5).
    %y = 2 + trnd(3, [t 1]);     % Simulate Student t (mue,dof=3) random numbers            

    m1(i) = sum(y.^1)/t;                                                      
    m2(i) = sum(y.^2)/t;                                                       
    m3(i) = sum(y.^3)/t;                                                        
    m4(i) = sum(y.^4)/t;                                                      

end
format ShortG

disp( [' Sample size    = ' num2str(t) ] );
disp( ' ' );

disp( ['Mean of m1     = ' num2str(mean(m1)) ] );
disp( ['Variance of m1 = ' num2str(std(m1)^2) ] );
disp( ' ' );
disp( ['Mean of m2     = ' num2str(mean(m2)) ] );
disp( ['Variance of m2 = ' num2str(std(m2)^2) ] );
disp( ' ' );
disp( ['Mean of m3     = ' num2str(mean(m3)) ] );
disp( ['Variance of m3 = ' num2str(std(m3)^2) ] );
disp( ' ' );
disp( ['Mean of m4     = ' num2str(mean(m4)) ] );
disp( ['Variance of m4 = ' num2str(std(m4)^2) ] );

