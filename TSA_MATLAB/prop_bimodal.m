%=========================================================================
%**
%**     Program to demonstrate multiple roots of the bivariate normal model
%**
%=========================================================================

clear all;
clc;


% This commented out code generates the data from scratch
% RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234) )

% t   = 4;                              % Sample size                                 
% rho = 0.6;
% u1  = randn(t,1);                     % Generate independent N(0,1)                 
% u2  = randn(t,1);
% x   = u1;                             % Generate dependent normal random numbers   
% y   = rho*u1 + sqrt(1 - rho^2)*u2;

% Take the GAUSS data for consistency
x =   [-0.60303846 ;
       -0.098331502; 
       -0.15897445 ;
       -0.65344600 ]; 

y =   [ 0.15367001; 
       -0.22971064; 
        0.66821992; 
       -0.44328369 ];
   
sxy = mean(x.*y);
sxx = mean(x.^2);
syy = mean(y.^2);
lnl = zeros(199,1);                   % log-likelihood                    
g   = zeros(199,1);                   % gradients                                 
rho = -0.99:0.01:0.99;                % grid values of rho      

tg = length( rho );

for i = 1:tg

    lnl(i) = -log(2*pi) - 0.5*log(1 - rho(i)^2) ...
             - 0.5*(sxx - 2*rho(i)*sxy + syy)/(1 - rho(i)^2);
    g(i)   = rho(i)*(1 - rho(i)^2) + (1 + rho(i)^2)*sxy - rho(i)*(sxx + syy);

end


%**************************************************************************
%**
%**     Generate graph
%**
%**************************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

subplot(1,2,1)
plot(rho,zeros(tg,1),'-k',rho,g,'-k','LineWidth',1);
title('(a) Gradient');
ylabel('$G_T(\rho)$');
xlabel('$\rho$');
axis([-1.0 1.0 -0.6 0.6])
set(gca,'YTick',[-0.6 -0.4 -0.2 0.0 0.2 0.4 0.6]);
set(gca,'XTick',[-1.0 -0.5 0.0 0.5 1.0]);
box off;

subplot(1,2,2)
plot(rho(5:end-5),lnl(5:end-5),'-k','LineWidth',1);
title('(b) Log-likelihood function');
ylabel('$L_T(\rho)$');
xlabel('$\rho$');
axis([-1.0 1.0 -3.0 -1.5])
set(gca,'YTick',[-3.0 -2.5 -2.0 -1.5]);

box off;

%laprint(1,'binormal','options','factory'); 