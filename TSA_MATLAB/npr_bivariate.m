%=========================================================================
%
%    Gaussian bivariate nonparametric kernel regression example 
%    
%=========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) );

% Simulate the nonlinear model
n = 2;          % Number of explanatory variables
t = 20000;       % Number of observations 

xt = randn(t,n);
ut = 0.1*randn(t,1);
yt = xt(:,1).*xt(:,2).^2 + ut;

% Set up a grid between -3 and 3
xmesh = (-3:0.2:3)';

% Need to work out all possible points at which to compute bivariate
% density on this mesh and express results as vectors
[x1m,x2m] = meshgrid(xmesh);

x = [reshape( x1m,[],1 ) reshape( x2m,[],1 )];

% Estimate density using product kernel
fac = 1/(4.0+n);
h   = std( xt )./(t.^fac);
ph  = prod( h );

ker  = zeros(n,1);
fx   = zeros(length(x),1);
fxy  = zeros(length(x),1);
pker = zeros(t,n);
for j = 1:length(x)
    
    for i = 1:t;
        
        for p = 1:n
            ker(p) = normpdf( (x(j,p) - xt(i,p))/h(p) );
        end 
        pker(i,1)  = prod( ker );
        pker(i,2) = prod( ker ).*yt(i);
    end
    fx(j)  = mean( pker(:,1) )/ph;
    fxy(j) = mean( pker(:,2) )/ph;
end

mx  = fxy./fx; 


%**************************************************************************
%**
%**     Generate graphs
%**
%**************************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

%--------------------------------------------------------%
% Panel (a)
subplot(1,2,1);
mesh( x1m,x2m,reshape( x(:,1).*(x(:,2).^2),length(x1m),length(x2m) ), ...
      'EdgeColor','black');
title('(a) True Surface');
xlabel('$x_{1,t}$');
ylabel('$x_{2,t}$');
zlabel('$f(x_{1,t},x_{2,t})$');
set(gca,'ztick',[])
axis tight
grid 'off'
box 'off'



%--------------------------------------------------------%
% Panel (b)
subplot(1,2,2);
mesh( x1m,x2m,reshape( mx,length(x1m),length(x2m) ),'EdgeColor','black');
title('(b) Estimated Surface');
xlabel('$x_{1,t}$');
ylabel('$x_{2,t}$');
zlabel('$f(x_{1,t},x_{2,t})$');
set(gca,'ztick',[])
axis tight
grid 'off'
box 'off'

%laprint(1,'bivariate','options','factory');

