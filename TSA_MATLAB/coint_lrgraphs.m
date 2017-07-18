%=========================================================================
%
%   Program to plot US macroeconomic data and 
%   to identify long-run properties of economic models
%=========================================================================

% Load perminent income data (total period is 1947q1 to 2010q2)
load permincome.mat

% Select desired sample (1984Q1 to 2005Q4)
c    = rcpc(150:237);
y    = rypc(150:237);

disp('OLS estimates of the permanent income equation (1984Q1 to 2005Q4)');
disp( [ ones(length(y),1) log(y)]\log(c)) ;


%**********************************************************************
%***
%***     Generate graph
%***
%**********************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;
    
dvec = seqa( 1984,1/4,length(c) );
%--------------------------------------------------------%
% Panel (a)
subplot(1,2,1)
plot( dvec,log(c),'-k',  ...
      dvec,log(y),'--k');
title('(a)')
ylabel('$lrc_t$, $lry_t$');
ylim([9.8 10.4] );
set(gca,'ytick',[ 9.8 10 10.2 10.4]);
box off


%--------------------------------------------------------%
% Panel (b)
subplot(1,2,2)
scatter( log(y),log(c),'.k');
title('(b)')
ylabel('$lrc_t$');
xlabel('$lry_t$');
box off
ylim([9.8 10.4] );
set(gca,'ytick',[9.8 10 10.2 10.4]);
xlim([9.8 10.4] );
set(gca,'xtick',[9.8 10 10.2 10.4]);

 
% Print the tex file to the relevant directory
%laprint(1,'permincome','options','factory');


clear all;

% Load money demand data (1959q1 to 2005q4)  
load moneydemand.mat


lrm    = log(m2./cpi);
lry    = log(gdp./cpi);
spread = tbill/100 - fedfunds/100;

x    = [ ones(length(lry),1) lry spread];
bhat = x\lrm;

disp('OLS estimates of the money demand (1959Q1 to 2005Q4)');
disp( bhat ) ;

u = lrm - x*bhat;


%**********************************************************************
%***
%***     Generate graph
%***
%**********************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(2);
clf;
    

%--------------------------------------------------------%
% Panel (a)
subplot(1,2,1)
scatter3(lry,spread,lrm,5,'k','filled');
title('(a)')
xlabel('$lry_t$');
ylabel('Spread');
zlabel('$lrm_t$');
box off

%--------------------------------------------------------%
% Panel (b)
subplot(1,2,2)
plot(seqa(1959,1/4,188),u,'-k');
title('(b)')
ylabel('Residual');
box off
ylim([-0.2 0.2]);
xlim([1959 2005]);
 
% Print the tex file to the relevant directory
%laprint(2,'moneydd','options','factory');


% Load interest rate data (March 1962 to September 2010)
load usmacro.mat
	

%**********************************************************************
%***
%***     Generate graph
%***
%**********************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(3);
clf;
    
dvec = seqa( 1962,1/4,length(r1yr) );
%--------------------------------------------------------%
% Panel (a)
subplot(1,2,1)
plot( dvec,r1yr,'-k',  ...
      dvec,r5yr,'--k', ...
      dvec,r10yr,':k');
title('(a)')
ylabel('Yield $\%$');
box off
axis tight


%--------------------------------------------------------%
% Panel (b)
subplot(1,2,2);
scatter3(r10yr,r5yr,r1yr,5,'k','filled');
title('(b)')
xlabel('10 year');
ylabel('5 year');
zlabel('1 year');
view([-69.5 10])
box off

 
% Print the tex file to the relevant directory
%laprint(3,'spread','options','factory');


