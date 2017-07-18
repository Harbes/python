%=========================================================================
%
%     Program to demonstrate the distribution of the t-statistic
%     from a Student t distribution with degrees of freedom of nu={1,2,3}.
%
%=========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',12) )

t = 500;                              
r = 5000;                            

%  Generate Student t (nu=1) ie Cauchy distribution     
nu     = 1;
chi1   = randn(t,r).^2;
rstud1 = randn(t,r)./sqrt( chi1./nu );

z1 = sqrt(t)*mean(rstud1)./std(rstud1);

figure(1)
histfit(z1,41);
title('Student t with nu=1')
xlabel('Midpoint')
ylabel('Frequency')

%  Generate Student t (nu=2) 
nu     = 2;
chi2   = randn(t,r).^2 + randn(t,r).^2;
rstud2 = randn(t,r)./sqrt( chi2./nu );

z2 = sqrt(t)*mean(rstud2)./std(rstud2);

figure(2)
histfit(z2,41);  
title('Student t with nu=2')
xlabel('Midpoint')
ylabel('Frequency')

%  Generate Student t (nu=3)  
nu     = 3;
chi3   = randn(t,r).^2 + randn(t,r).^2 + randn(t,r).^2;
rstud3 = randn(t,r)./sqrt( chi3./nu );

z3 = sqrt(t)*mean(rstud3)./std(rstud3);

figure(3)
histfit(z3,41);                         
title('Student t with nu=3')
xlabel('Midpoint')
ylabel('Frequency')
