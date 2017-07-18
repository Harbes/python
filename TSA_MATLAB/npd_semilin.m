%=========================================================================
%
%   Semi-parametric Look Ahead Estimator (LAE): 
%   linear autoregressive model
%
%=========================================================================

clear all
clc

RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) )

% Parameters
mu    = 1.0;
phi   = 0.5;
sig2  = 1.0;
t     = 5;

% Grid
y = -4:0.1:6;
n = length(y);

%  Simulate AR(1) model   
yt = mu*ones(t,1);

for i = 2:t

    yt(i) = mu + phi*(yt(i-1)) + sqrt(sig2)*randn(1,1);


end

% Generate the stationary distribution using the LAE 
m = mu + phi*yt;
s = sqrt(sig2);
    
f = zeros(length(yt),length(y));
for j = 1:length(yt)
    
   for i = 1:length(y)    
        
       z    = (y(i) - m(j))/s;
       f(j,i)= (1/sqrt(2*pi*s^2))*exp(-0.5*(z^2));

   end
end

f_lae = mean( f );

% Generate the true stationary distribution using the analytical expression      
m = mu/(1 - phi);
s = sqrt( sig2/(1 - phi^2));
z = ( y - m )/ s ;

f_true = (1/sqrt(2*pi*s^2)).*exp(-0.5*(z.^2));  


% Kernel density estimate of the stationary distribution    
h     = 1.06*std(yt)*t^(-1/5);  
    
fx = zeros(length(y),1);
for i = 1:length(y);
    
    z      = ((y(i) - yt)./h);    
    fx(i)  = mean( normpdf(z)./h );

end



%	f_lae_components = ( (1/sqrt(2*pi*s^2)).*exp(-0.5*(z.^2)) )';  

%**********************************************************************
%***
%***     Generate graph
%***
%**********************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

plot(y,f_lae,'--k','LineWidth',1.0);
hold on;
plot(y,f_true,'-k','LineWidth',0.5);
plot(y,fx,':k');
title('True (solid), LAE (dashed) and Kernel Densities');
ylabel('$f(y)$');
xlabel('$y$');
box off;





