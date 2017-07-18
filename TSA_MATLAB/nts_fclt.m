%=========================================================================
%
%   Program to demonstrate the properties of the 
%   Functional Central Limit Theorem: 
%   Standardization of a random walk
%
%=========================================================================
clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) )
       
ndraws = 50000;
t      = 500;

alpha = 0.0;
phi   = 1.0;
sig2  = 1.0;
y0    = 0.0;

yts1 = zeros(ndraws,1);    % s = 1/4
yts2 = zeros(ndraws,1);    % s = 1/2
yts3 = zeros(ndraws,1);    % s = 1

for i = 1:ndraws

    % Random walk with first observation discarded
    y  = trimr( recserar(alpha + sqrt(sig2)*randn(t+1,1),y0,phi),1,0 );   

    
    yts1(i) = y(round(t*0.25))*t^(-0.5)/sqrt(sig2);                                                
    yts2(i) = y(round(t*0.5))*t^(-0.5)/sqrt(sig2);                                              
    yts3(i) = y(end)*t^(-0.5)/sqrt(sig2);                                                             

end



disp(['Sample size = ',num2str(t)] );
disp(' ');
disp(['Sample mean of yts1 (T/4)               = ', num2str( mean(yts1) )] );
disp(['Theoretical mean of yts1 (T/4)          = ', num2str( 0.0)] );
disp(['Sample variance of yts1                 = ', num2str( std(yts1)^2)] );
disp(['Theoretical variance of yts1            = ', num2str( 1/4)] );

disp(' ');

disp(['Sample mean of yts2 (T/2)               = ', num2str( mean(yts2))] );
disp(['Theoretical mean of yts2 (T/2)          = ', num2str( 0.0)] );
disp(['Sample variance of yts2                 = ', num2str( std(yts2)^2)] );
disp(['Theoretical variance of yts2            = ', num2str( 1/2)] );

disp(' ');

disp(['Sample mean of yts3 (T)                 = ', num2str( mean(yts3))] );
disp(['Theoretical mean of yts3 (T)            = ', num2str( 0.0)] );
disp(['Sample variance of yts3                 = ', num2str( std(yts3)^2)] );
disp(['Theoretical variance of yts3            = ', num2str( 1/1)] );
disp(' ');


xi = seqa(-5,0.1,101);

fhat1 = ksdensity(yts1,xi);
fhat2 = ksdensity(yts2,xi);
fhat3 = ksdensity(yts3,xi);
fnorm = normpdf( xi );


%**********************************************************************
%***
%***     Generate graph
%***
%**********************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;
    
%--------------------------------------------------------%
plot( xi,fnorm,'-k',  ...
      xi,fhat1,'--k', ...
      xi,fhat2,'-.k', ...
      xi,fhat3,':k');
box off
xlim( [-4 4]);
set(gca,'xtick',[-4 -3 -2 -1 0 1 2 3 4]);
xlabel('$Y_{T}(s)$');
ylabel('$f(Y_{T}(s))$');
 
% Print the tex file to the relevant directory
%laprint(1,'fclt','options','factory');

