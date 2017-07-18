%=========================================================================
%
%   Semi-parametric Look Ahead Estimator (LAE): 
%   nonlinear threshold autoregressive model
%
%=========================================================================

clear all
clc

RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) )

% Parameters
theta = 0.6;
t     = 5;


%  Simulate tar(1) model   
yt = theta*ones(t,1);

for i = 2:t

    yt(i) = theta*abs(yt(i-1)) + sqrt(1 - theta^2)*randn(1,1);


end

% or set values to those in text

yt = [ 0.600; -0.279; -1.334; 0.847; 0.894 ];

% Grid
y = -4:0.1:4;
n = length(y);

% Generate the stationary distribution using the LAE 
s = sqrt(1 - theta^2);
    
f = zeros(length(yt),length(y));
for j = 1:length(yt)
    
   for i = 1:length(y)    
        
       z    = (y(i) - theta*yt(j))/s;
       f(j,i)= (1/sqrt(2*pi*s^2))*exp(-0.5*(z^2));

   end
end

f_lae = mean( f );



% Generate the true stationary distribution using the analytical expression      
delta  = theta/sqrt(1 - theta^2);
f_true = 2*normpdf(y).*normcdf(delta*y);    

% % Kernel density estimate of the stationary distribution    **/  
h     = 1.06*std(yt)*t^(-1/5);  
    
fx = zeros(length(y),1);
for i = 1:length(y);
    
    z      = ((y(i) - yt)./h);    
    fx(i)  = mean( normpdf(z)./h );

end

%**********************************************************************
%***
%***     Generate graphs
%***
%**********************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

%--------------------------------------------------------%
% Panel (a)
subplot(1,2,1);
plot(y,f(1,:),'--k','LineWidth',1.0);
hold on;
plot(y,f(3,:),'-k','LineWidth',0.5);
plot(y,f(5,:),':k');
hold off;
title('(a) Component Transitional Densities');
ylabel('$f(y)$');
xlabel('$y$');
% legend('Estimated','True','Location','NorthWest')
% legend boxoff;
box off;
hold off;

%--------------------------------------------------------%
% Panel (b)
subplot(1,2,2);
plot(y,f_lae,'--k','LineWidth',1.0);
hold on;
plot(y,f_true,'-k','LineWidth',0.5);
plot(y,fx,':k');
title('(b) True, LAE and Kernel');
ylabel('$f(y)$');
xlabel('$y$');
% legend('Estimated','True','Location','NorthWest')
% legend boxoff
box off;

%laprint(1,'seminonlin','options','factory');



