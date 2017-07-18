%=========================================================================
%
%     Simulation of nonlinear and linear exponential models.
%
%=========================================================================

clear all
clc

RandStream.setDefaultStream( RandStream('mt19937ar','seed',1) )

% Simulate data 
t   = 50;        		     
b0  = 1.0;
b1  = 0.05;
sig = 0.5;
u   = sig*randn(t,1);
x   = (1:1:t)';
y1  = b0*exp( b1*x + u);
y2  = b0*exp( b1*x ) + u;



%***************************************************************************
%***
%***     Generate graphs
%***
%***************************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

%--------------------------------------------------------%
% Panel (a)
subplot(1,2,1);
plot(x,y1,'-k',...
     x,y2,'-.k',...
     'LineWidth',0.75);

title('(a) Levels');
xlabel('$x_t$');
ylabel('$y_t$');
set(gca,'XTick',0:10:50);
set(gca,'YTick',0:5:20);
xlim([0,50]);
ylim([0,20]);
%set(gca,'LineWidth',1);
%legend('$y_{1,t}$','$y_{2,t}$','Location','NorthWest');
box off;
%legend('boxoff');

%--------------------------------------------------------%
% Panel (b)
subplot(1,2,2);
plot(x,log(y1),'-k',...
     x,log(y2),'-.k',...
     'LineWidth',0.75);

title('(b) Logs');
xlabel('$x_t$');
ylabel('$\log y_t$');
set(gca,'XTick',0:10:50);
set(gca,'YTick',-1:1:4);
xlim([0,50]);
ylim([-1,4]);
%set(gca,'LineWidth',1);
%legend('$\log y_{1,t}$','$\log y_{2,t}$','Location','NorthWest');
box off;
legend('boxoff');

 
%laprint(1,'nonlinexpsim','options','factory');


