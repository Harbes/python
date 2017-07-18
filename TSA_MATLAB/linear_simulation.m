%=========================================================================
%
%    Generate simulated data by simulating the reduced form.
%    The set of equations is given by the bivariate system
%        y1t = beta1*y2t + alpha1*x1t + u1t
%        y2t = beta2*y1t + alpha2*x2t + u2t
%     where E[ut'ut] = omega
%
%=========================================================================


clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',123457) )

t = 500; 

beta1  = 0.6; 
beta2  = 0.2; 
alpha1 = 0.4;
alpha2 = -0.5;

omega = [1 0.5; 0.5 1.0];

% Construct population parameter matrices     
B = [1 -beta2; -beta1 1]; 
A = [-alpha1 0; 0 -alpha2];

% Construct exogenous variables                       
x = [10*randn(t,1) 3*randn(t,1)];

% Construct structural disturbances                   
u = randn(t,2)*chol(omega);

% Construct reduced form parameters                   
invB = inv(B);
phi  = -A*invB;

disp('Inverse of b');
disp(invB);
disp('Reduced form parameter matrix');
disp(phi);

% Construct reduced form disturbances                 
v = u*invB;

% Simulate the model by simulating the reduced form   
y = zeros(t,2);

for i=1:t
    y(i,:) = -x(i,:)*A*invB + v(i,:);
end


% ***     Generate graphs  ***

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

t = 1:1:500;

%--------------------------------------------------------%
% Panel (a)
subplot(2,2,1);
plot(t,y(:,1),'-k');
title('(a)');
xlabel('$t$');
ylabel('$y_{1,t}$');
set(gca,'XTick',0:100:500);
set(gca,'YTick',-10:5:10);
xlim([0,500]);
ylim([-10,10]);
box off;
%set(gca,'LineWidth',1);

%--------------------------------------------------------%
% Panel (b)
subplot(2,2,2);
plot(t,y(:,2),'-k');
title('(b)');
xlabel('$x_t$');
ylabel('$y_{2,t}$');
set(gca,'XTick',0:100:500);
set(gca,'YTick',-10:5:10);
xlim([0,500]);
ylim([-10,10]);
box off;
%set(gca,'LineWidth',1);

%--------------------------------------------------------%
% Panel (c)
subplot(2,2,3);
plot3(x(:,1),y(:,1),y(:,2),'.k','MarkerSize',5);
title('(c)');
xlabel('$x_{1,t}$');
ylabel('$y_{1,t}$');
zlabel('$y_{2,t}$');
set(gca,'XTick',-10:5:10);
set(gca,'YTick',-10:5:10);
xlim([-10,10]);
ylim([-10,10]);
%set(gca,'LineWidth',1);

%--------------------------------------------------------%
% Panel (d)
subplot(2,2,4);
plot3(x(:,2),y(:,2),y(:,1),'.k','MarkerSize',5);
title('(d)');
xlabel('$x_{2,t}$');
ylabel('$y_{2,t}$');
zlabel('$y_{1,t}$');
set(gca,'XTick',-10:5:10);
set(gca,'YTick',-10:5:10);
xlim([-10,10]);
ylim([-10,10]);


%laprint(1,'syssim','options','factory');
