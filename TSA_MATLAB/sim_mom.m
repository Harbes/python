% ========================================================================
%
%     Method of moments estimation of a first order MA model. 
%
% ========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) );

t     = 250;        % Sample size             
theta = 0.5;        % MA1 Population parameter

% Generate the data for a MA(1) model      
u  = randn(t,1);
y  = u(2:end) - theta*u(1:end-1);

% Estimate the first order autocorrelation coefficient
y   = y - mean(y);
rho = y(2:end)\y(1:end-1);


% Estimate theta using the method of moments estimator   
b_mom = ( -1 + sqrt(1 - 4*rho^2) ) / (2*rho);


disp(' ')
disp(['Sample size                        = ', num2str(t) ]);
disp(['True population parameter (theta)  = ', num2str(theta) ]);
disp(['Method of moment estimate          = ', num2str(b_mom) ]);
disp(['True population parameter (AR(1))  = ', num2str(-theta/(1+theta^2)) ]);
disp(['First order AR(1)                  = ', num2str(rho) ]);


%**************************************************************************
%**
%**     Generate graph
%**
%**************************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;


plot(1:1:249,y,'-k');
xlabel('$t$');
ylabel('$y_t$');
axis tight;
box off;
set(gca,'XTick',0:50:250);
set(gca,'YTick',-4:2:4);
xlim([0,250]);
ylim([-4,4]);

%laprint(1,'ma1','options','factory');



