%=========================================================================
%
%   Program to demonstrate the distribution of the 
%   estimates of a gamma regression model
%
%=========================================================================
clear all
clc
    
ndraws = 5000;                                                         

% Parameters
b0    = [ 1  2 ];                        
rho   = 0.25;                           
alpha = 0.1;

% For the sample size of T = 10     
RandStream.setDefaultStream( RandStream('mt19937ar','seed',1) );

t  = 10;                           
x  = [ones(t,1)  randn(t,1)];
z1 = zeros(ndraws,1);
z2 = zeros(ndraws,1);

for i = 1:ndraws
        
    u  = alpha*gamrnd(rho,1,[t 1]) - rho*alpha;                      
    y  = x*b0' + u;
    b  = x\y; 
    e  = y - x*b;
    s2 = e'*e/t;
    vc = s2*inv(x'*x);
   
    z1(i) = (b(1) - b0(1))/sqrt(vc(1,1));
    z2(i) = (b(2) - b0(2))/sqrt(vc(2,2));

end

% For the sample size of T = 100 
RandStream.setDefaultStream( RandStream('mt19937ar','seed',1) );

t  = 100;                              
x  = [ones(t,1)  randn(t,1)];
z3 = zeros(ndraws,1);
z4 = zeros(ndraws,1);

for i = 1:ndraws
        
    u  = alpha*gamrnd(rho,1,[t 1]) - rho*alpha;                      
    y  = x*b0' + u;
    b  = x\y; 
    e  = y - x*b;
    s2 = e'*e/t;
    vc = s2*inv(x'*x);
   
    z3(i) = (b(1) - b0(1))/sqrt(vc(1,1));
    z4(i) = (b(2) - b0(2))/sqrt(vc(2,2));

    end

% Plot the results
bin = -5.25:0.5:4.75;
    
figure(1)

subplot(2,2,1)  
hist(z1,bin); 
title('(a) T=10, $\beta_0$')
axis([-5 5 0 1200])


subplot(2,2,2)
hist(z2,bin);                         
title('(a) T=10, $\beta_1$')

subplot(2,2,3)  
hist(z3,bin);                         
title('(a) T=100, $\beta_0$')

subplot(2,2,4)
hist(z4,bin);                         
title('(a) T=100, $\beta_1$')

