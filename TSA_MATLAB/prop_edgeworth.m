%=========================================================================
%
%   Program to compare the edgeworth expansion of the finite sample distribution
%   and the asymptotic distribution (exponential distribution) 
%
%=========================================================================

clear all
clc

t = 5;           
s = -2:1:2;   

% Finite sample value (based on a gamma distribution)
f_finite = 1 - gamcdf( t./(1 + s/sqrt(t)),t,1 );             

% Asymptotic value (based on the normal distribution)
f_asy = normcdf( s );                

% Compute Edgeworth expansion values 
h2 = s.^2 - 1;               
h3 = s.^3 - 3*s;
h5 = s.^5 - 10*s.^3 + 15*s;

f_ew = normcdf(s) - normpdf(s).*( (1 + (2/3)*h2)/sqrt(t) + ( (5/2)*s + (11/12)*h3 + (2/9)*h5 )/t );   

% Print results
format ShortG
disp( [' Sample size    = ' num2str(t) ] );
disp( ' ' );
disp( '       s       Finite    Edgeworth  Asymptotic ');
disp([ s' f_finite' f_ew' f_asy'] );

 