%=========================================================================
%
%   Comparing the true likelihood (lnl0) and the quasi-likelihood (lnl) 
%
%=========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',12) )
t = 10;

% Exponential versus normal
% -------------------------
mu = 1;
y  = -mu.*log(1 - rand(t,1));                 
y = round(y*1000)/1000;            % Round y to three decimal places



% True likelihood (exponential)
theta0 = seqa(0.01,0.1,521);
lnl0   = -log(theta0) - mean(y)./theta0;         

% Quasi-likelihood (normal with unit variance)
theta = seqa(-5,0.1,201);
sig2  = 1;
lnl   = zeros(length(theta),1);

for i = 1:length(theta)
    
    lnl(i)   = -0.5*log(2*pi) - 0.5*log(sig2) ...
               - 0.5*mean( (y - theta(i) ).^2 )./sig2;

end



%********************************************************************
%***
%***     Generate graph
%***
%********************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

subplot(1,2,1)
plot(theta0,lnl0,'-k')
ylim([-2 -1])
xlim([-1  3])
title( '(a) Exponential - Normal')
ylabel( ' $\ln L_T(\theta)$ ') 
xlabel( ' $\theta$ ') 
hold on
plot(theta,lnl,'--k')
hold off
   
% Negative Binomial versus Poisson
% --------------------------------
mu = 5;
p  = 0.5;
y = nbinrnd(mu,p,[t,1]);                                         


% True likelihood (negative binomial)
theta0 = seqa(0.01,0.1,201);
lnl0   = zeros( length(theta0),1 );

for i = 1:length(theta0)
    
    lnl0(i) = mean(log(gamma(y + theta0(i)))) - mean(log(gamma(y + 1))) ...
            - log(gamma(theta0(i))) + theta0(i)*log(1-p) + mean(y)*log(p) ;         

end

% Quasi-likelihood (Poisson)

theta = seqa(0.01,0.1,201);
lnl = mean(y)*log(theta) - theta - mean(log(factorial(y)));  
  
subplot(1,2,2)
plot(theta0,lnl0,'-k')
ylim([-6 -2])
xlim([ 0  10]) 
title( '(b) Negative Binomial - Poisson')
ylabel( ' $\ln L_T(\theta)$ ') 
xlabel( ' $\theta$ ') 
hold on
plot(theta,lnl,'--k')
hold off


% Print the tex file to the relevant directory
%laprint(1,'qmlegraphs','options','factory');

% Student t versus Normal
% -----------------------
mu  = 5;
sig = 1;
v   = 10;
t   = 10;

% Standardized Student t random numbers
y = mu + sig*sqrt((v-2)/v)*tinv(rand(t,1),v);                    


% True likelihood (Student t)
theta0 = seqa(-5,0.1,201);
z      = zeros(t,length(theta0));
for i = 1:length(theta0)
    
    z(:,i) = y-theta0(i);
end

const  = gamma( (v+1)/2 ) / ( sqrt(pi*(v-2)) * gamma( v/2 ) );               

lnl0 = log(const) - 0.5*log(sig^2) - 0.5*(v+1)*mean( log( 1 + (z.^2)/(v-2) ) );            

% Quasi-likelihood (normal)
theta = seqa(-5,0.1,201);
for i = 1:length(theta0)    
    z(:,i) = y-theta0(i);
end
lnl = -0.5*log(2*pi) - 0.5*log(sig^2) - 0.5*mean( z.^2 );         

figure(2)
plot(theta0,lnl0,'-k')
title( '(b) t - Normal')
ylabel( ' $\ln L_T(\theta)$ ') 
xlabel( ' $\theta$ ') 
hold on
plot(theta,lnl,'--k')
hold off

 
% Poisson versus normal
%----------------------
mu = 5;
t = 10;

% Poisson random numbers 
y = poissrnd(mu,t,1);                                         

% True likelihood (Poisson)
theta0 = seqa(0.01,0.1,201);

lnl0 = mean(y)*log(theta0) - theta0 - mean(log(factorial(y)));  

% Quasi-likelihood (normal with variance = mu)
theta = seqa(-5,0.1,201);
sig2  = mu; 
for i = 1:length(theta0)    
    z(:,i) = y-theta0(i)/sig2;
end

lnl = -0.5*log(2*pi) - 0.5*log(sig2) - 0.5*mean( z );


figure(3)
plot(theta0,lnl0,'-k')
title( ' Poisson - Normal')
ylabel( ' $\ln L_T(\theta)$ ') 
xlabel( ' $\theta$ ') 
hold on
plot(theta,lnl,'--k')
hold off

