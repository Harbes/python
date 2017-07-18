%=========================================================================
%
%     Semiparametric kernel regression  
%
%=========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) );

% Simulate the model
t   = 500;
ut  = 0.1*randn(t,1);                       % N(0,0.1^2) 

x3t = sortrows( -2 + 4*rand(t,1) );         % U(-2,2) 
x1t = 0.5*x3t + randn(t,1);                 % N(0.5*x3t,1)  
x2t = trnd( 4,t,1 );                        % Standardised Student t (0,1,4)  

f   = 0.3*exp( -4*(x3t+1).^2 ) + 0.7*exp( -16*(x3t-1).^2 );
yt  = 2*x1t + 1*x2t + f + ut;

% Kernel regression of yt on x3t 
h  = 0.05;
fx = zeros(t,1);
fxy = zeros(t,1);
for i = 1:t;
    z      = ((x3t(i) - x3t)./h);    
    fx(i)  = mean( normpdf(z)./h );
    fxy(i) = mean( normpdf(z).*yt./h );
end
mxy = fxy ./ fx;

% Kernel regression of x1t on x3t 
for i = 1:t;
    z      = ((x3t(i) - x3t)./h);    
    fx(i)  = mean( normpdf(z)./h );
    fxy(i) = mean( normpdf(z).*x1t./h );
end
mx1 = fxy./fx;

% Kernel regression of x2t on x3t
for i = 1:t;
    z      = ((x3t(i) - x3t)./h);    
    fx(i)  = mean( normpdf(z)./h );
    fxy(i) = mean( normpdf(z).*x2t./h );
end
mx2 = fxy./fx;

% Now run linear regression
y = yt - mxy;
x = [(x1t - mx1)  (x2t - mx2)];
b = inv(x'*x)*x'*y;

disp('Parameter estimates from semi-parametric model');
disp('**********************************************');
disp(b');


% Reconstruct f 

fhat = mxy - b(1)*mx1 - b(2)*mx2;


% OLS regression of yt on x1t, x2t and x3t 
y = yt;
x = [x1t  x2t  x3t];
b = x\y;

disp('Parameter estimates from regressing yt on x1t, x2t and x3t');
disp('**********************************************************');
disp(b');


% OLS regression of yt on x1t, x2t and x3t^2 
y = yt;
x = [x1t  x2t  (x3t.^2)];
b = x\y;

disp('Parameter estimates from regressing yt on x1t, x2t and x3t^2');
disp('************************************************************');
disp(b');


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
plot(1:1:t,yt,'-k');
title('(a) Data');
xlabel('$t$');
ylabel('$y_t$');
axis tight
box 'off'

%--------------------------------------------------------%
% Panel (b)
subplot(1,2,2);
plot(x3t,f,'-.k',...
     x3t,fhat,'-k');
title('(b) Estimated Functional Form');
xlabel('$f(x_{3,t})$');
ylabel('$x_{3,t}$');
axis tight
box 'off'

%laprint(1,'semiparam','options','factory');

