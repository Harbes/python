%**************************************************************************
%**
%**     Program to generate Figure 1a of Stock, Wright and Yogo 
%**     weak instrument paper (JBES, 2002)
%**
%**************************************************************************
clear all;
clc;

state = 1234;
rand('twister', state);
randn('state', state);


N     = 10000;          % Number of replications  
T     = 5;              % Sample size            
beta  = 0.0;            % Parameter values      
phi   = 0.25;
sig11 = 1.0;
sig22 = 1.0;
sig12 = 0.99;
omega = [sig11 sig12 ; sig12 sig22];
rho   = sig12/sqrt(sig11*sig22);

%**************************************************************************
%**
%**     Generate instruments (scaled) and do Monte Carlo replications
%**
%**************************************************************************

x = rand(T,1);                 
%x = x/sqrt(x'*x);

biv = zeros(N,1);

for i = 1:N

    %u = mvnrnd([0 0],omega,T);
    u  = randn(T,length(omega))*chol(omega);
    w = x*phi  + u(:,1);                                          % Reduced form equation   
    y = w*beta + u(:,2);                                          % Structural equation     

    tmp    = inv(x'*x);
    tmp1   = x*tmp*x';
    tmp2   = w'*tmp1;
    biv(i) = inv(tmp2*w)*(tmp2*y);      % IV estimates            
end

%**************************************************************************
%**
%**     Generate graph
%**
%**************************************************************************



xi = -2.0:0.001:2.0;
f  = ksdensity(biv,xi);

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

plot(xi,f,'-k');
xlabel('$\beta_{\text{IV}}$');
ylabel('$f(\beta_{\text{IV}})$');
box off;

%laprint(1,'weak','options','factory');
 