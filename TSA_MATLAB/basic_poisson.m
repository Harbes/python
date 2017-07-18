%=========================================================================
%
%   Program to estimate an Poisson model and plot the 
%   log-likelihood function.
%
%=========================================================================

clear all;
clc;

% Data
y = [ 8, 3, 4 ];        
% y = [ 6, 2, 3, 1 ] ;    % Data used in Poisson exercies

t    = length(y);      
theta = mean(y);
lnl_t = y*log(theta) - theta - log(factorial(y));
g_t   = y/theta - 1;
h_t   = -y/theta^2;

disp( ['Sum of y   = ', num2str( sum(y) )] );
disp( ['Mean of y  = ', num2str( theta )] );
disp(['Log-likelihood function = ', num2str(mean(lnl_t)) ] );

disp('');
disp( '       yt      lnlt       gt       ht  ');
disp( [y'   lnl_t'  g_t'  h_t' ]);


% ***    Generate graph   ***

theta = 0.01:0.01:15;
lnl   = zeros(length(theta),1);

for i = 1:length(theta)
    
    lnl(i)   = mean( y*log(theta(i)) - theta(i)  - log(factorial(y)) );
    
end

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

figure1 = figure;

% Create subplot
subplot1 = subplot(1,2,1,'Parent',figure1);
hold(subplot1,'all');
plot(theta,lnl,'Parent',subplot1,'Color',[0 0 0]);
title('(a) Log-likelihood function');
ylabel('$\ln L_T(\theta)$');
xlabel('$\theta$');
box off;


subplot2 = subplot(1,2,2,'Parent',figure1);
view(subplot2,[0 -90]);
hold(subplot2,'all');
bar(y,lnl_t,'Parent',subplot2);
title('(b) Log-density function');
ylabel('$\ln f(y_t;5)$','VerticalAlignment','bottom','Rotation',90,...
    'HorizontalAlignment','center');
xlabel('$y_t$','VerticalAlignment','cap','HorizontalAlignment','center');
set(gca,'xtick',[1 2 3 4 5 6 7 8 9 10])

laprint(figure1,'poisson','options','factory');
