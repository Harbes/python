%**************************************************************************
%**
%**     Program to generate graphs in the hypothesis testing chapter
%**
%**************************************************************************
clear all;
clc;

state = 123457;
rand('state', state);
randn('state', state);


%**************************************************************************
%**     
%**    Likelihood ratio statistic graph        
%**                                         
%**************************************************************************
T = 100;

theta0 = 0.0;                       % Population mean (unknown)       
sig2   = 1.0;                       % Population variance (known)     

y = theta0 + sqrt(sig2)*randn(T,1); % Simulate the data              

thetahat = mean(y);                 %    MLE: sample mean               

theta = (-0.02:0.001:0.02)';

tmp   = repmat( y,1,length(theta) );
tmp1  = repmat( theta',length(y),1 );
tmp3  = tmp - tmp1;
lf    = -0.5*T*log(2*pi) - 0.5*T*log(sig2) - 0.5*sum( tmp3.^2/sig2 );      %    Log likelihood                          
lfhat = -0.5*T*log(2*pi) - 0.5*T*log(sig2) - 0.5*sum( ( (y - thetahat').^2 )/sig2 );  %     Log-likelihood evaluated at thetahat    
lf0   = -0.5*T*log(2*pi) - 0.5*T*log(sig2) - 0.5*sum( ( (y - theta0').^2 )/sig2 );    %    Log-likelihood evaluated at theta0     

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;
figure(1)
plot(theta,lf,'-k','LineWidth',1.5)
set(gca,'YTick',[-130.4866 -130.4828])
set(gca,'YTickLabel',{'$\ln L(\theta_0)$','$\ln L(\widehat{\theta})$'})
set(gca,'XTick',[-0.01 -0.001])
set(gca,'XTickLabel',{'$\theta_0$','$\widehat{\theta}$'})
box off;


% Create line for theta hat
line([-0.001 -0.001],[-130.51 -130.4828],'Color','k','LineStyle',':');
line([-0.02  -0.001],[-130.4828 -130.4828],'Color','k','LineStyle',':');

% Create line for theta null
line([-0.01 -0.01],[-130.51 -130.4866],'Color','k','LineStyle',':');
line([-0.02  -0.01],[-130.4866 -130.4866],'Color','k','LineStyle',':');


%laprint(1,'C:\Stan\Current\Vancebook\Graphics\m_fig-testlr','options','factory'); 

%**************************************************************************
%**     
%**    Lagrange multiplier statistic graph        
%**                                         
%**************************************************************************



g   = 0.5*sum( tmp3/sig2 );                    % Gradient function 
g_thetahat = 0.5*sum( (y - thetahat')/sig2 );  % Gradient evaluated at thetahat      **/
g_theta0   = 0.5*sum( (y - theta0')/sig2 );    % Gradient evaluated at thetahat      **/

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(2);
clf;

plot(theta,g','-k','LineWidth',1.5);
hold on;
plot(theta,ones(length(g),1)*(-0.0166),':k');

% ylabel('$G(\theta)$');
% xlabel('$\theta$');
box off;
axis tight;

set(gca,'YTick',[-0.0166 0.4334])
set(gca,'YTickLabel',{'$G(\widehat{\theta})=0$','$G(\theta_0)$'})

set(gca,'XTick',[-0.01 -0.001])
set(gca,'XTickLabel',{'$\theta_0$','$\widehat{\theta}$'})

% Create line for theta hat
line([-0.001 -0.001],[-1.0666 -0.0166],'Color','k','LineStyle',':');

% Create line for theta null
line([-0.01 -0.01],[-1.0666 0.4334],'Color','k','LineStyle',':');
line([-0.02 -0.01],[0.4334 0.4334],'Color','k','LineStyle',':');

%laprint(2,'C:\Stan\Current\Vancebook\Graphics\m_fig-testlm','options','factory'); 
