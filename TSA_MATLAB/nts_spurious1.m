%=========================================================================
%
%   Program to demonstrate the spurious regression problem 
%   using correlation coefficients
%
%=========================================================================
clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',15) )

t = 100;                            
r = 10000;                          



% I(0) versus I(0)      
rca = zeros(r,1);

for i = 1:r
    
    y1   = trimr(recserar(randn(t+100,1),0.0,0.0),100,0);
    y2   = trimr(recserar(randn(t+100,1),0.0,0.0),100,0);
    rmat = corrcoef([ y1 y2 ]);
    rca(i) = rmat(1,2);
     
end

y = quantile( rca,[.01 .025 .05 .10 .50 .90 .95 .975 .99] );

disp( ' ' );
disp( 'The empirical distribution (stationary case)' );
disp( '********************************************' );
disp( [' 1.0 per cent      =',num2str(y(1)) ] );
disp( [' 2.5 per cent      =',num2str(y(2)) ] );
disp( [' 5.0 per cent      =',num2str(y(3)) ] );
disp( ['10.0 per cent      =',num2str(y(4)) ] );
disp( ['50.0 per cent      =',num2str(y(5)) ] );
disp( ['90.0 per cent      =',num2str(y(6)) ] );
disp( ['95.0 per cent      =',num2str(y(7)) ] );
disp( ['97.5 per cent      =',num2str(y(8)) ] );
disp( ['99.0 per cent      =',num2str(y(9)) ] );

% p   = [0.025 0.975 ];                      
% tmp = norminv(p,0,1)./sqrt(t);                        
% disp( ['95% confidence interval  ',num2str(tmp) ] );
    
 
% I(1) versus I(1)       
rcb = zeros(r,1);

for i = 1:r
    
    y1   = trimr(recserar(randn(t+100,1),0.0,1.0),100,0);
    y2   = trimr(recserar(randn(t+100,1),0.0,1.0),100,0);
    rmat = corrcoef([ y1 y2 ]);
    rcb(i) = rmat(1,2);
     
end

y = quantile( rcb,[.01 .025 .05 .10 .50 .90 .95 .975 .99] );

disp( ' ' );
disp( 'The empirical distribution (case 2)' );
disp( '********************************************' );
disp( [' 1.0 per cent      =',num2str(y(1)) ] );
disp( [' 2.5 per cent      =',num2str(y(2)) ] );
disp( [' 5.0 per cent      =',num2str(y(3)) ] );
disp( ['10.0 per cent      =',num2str(y(4)) ] );
disp( ['50.0 per cent      =',num2str(y(5)) ] );
disp( ['90.0 per cent      =',num2str(y(6)) ] );
disp( ['95.0 per cent      =',num2str(y(7)) ] );
disp( ['97.5 per cent      =',num2str(y(8)) ] );
disp( ['99.0 per cent      =',num2str(y(9)) ] );

% p   = [0.025 0.975 ];                      
% tmp = norminv(p,0,1)./sqrt(t);                       
% disp( ['95% confidence interval',num2str(tmp) ] );


% I(1) versus I(2)             
rcc = zeros(r,1);

for i = 1:r
    
    y1   = trimr(recserar(randn(t+100,1),0.0,1.0),100,0);
    y2   = trimr(recserar(randn(t+100,1),[0.0;0.0],[2.0;-1.0]),100,0);
    rmat = corrcoef([ y1 y2 ]);
    rcc(i) = rmat(1,2);
     
end

y = quantile( rcc,[.01 .025 .05 .10 .50 .90 .95 .975 .99] );

disp( ' ' );
disp( 'The empirical distribution (case 3)' );
disp( '********************************************' );
disp( [' 1.0 per cent      =',num2str(y(1)) ] );
disp( [' 2.5 per cent      =',num2str(y(2)) ] );
disp( [' 5.0 per cent      =',num2str(y(3)) ] );
disp( ['10.0 per cent      =',num2str(y(4)) ] );
disp( ['50.0 per cent      =',num2str(y(5)) ] );
disp( ['90.0 per cent      =',num2str(y(6)) ] );
disp( ['95.0 per cent      =',num2str(y(7)) ] );
disp( ['97.5 per cent      =',num2str(y(8)) ] );
disp( ['99.0 per cent      =',num2str(y(9)) ] );

% p   = [0.025 0.975 ];                      
% tmp = norminv(p,0,1)./sqrt(t);                       
% disp( ['95% confidence interval',num2str(tmp) ] );
 
% I(2) versus I(2)       
rcd = zeros(r,1);

for i = 1:r
    
    y1   = trimr(recserar(randn(t+100,1),[0.0;0.0],[2.0;-1.0]),100,0);
    y2   = trimr(recserar(randn(t+100,1),[0.0;0.0],[2.0;-1.0]),100,0);
    rmat = corrcoef([ y1 y2 ]);
    rcd(i) = rmat(1,2);
     
end

y = quantile( rcd,[.01 .025 .05 .10 .50 .90 .95 .975 .99] );

disp( ' ' );
disp( 'The empirical distribution (case 4)' );
disp( '********************************************' );
disp( [' 1.0 per cent      =',num2str(y(1)) ] );
disp( [' 2.5 per cent      =',num2str(y(2)) ] );
disp( [' 5.0 per cent      =',num2str(y(3)) ] );
disp( ['10.0 per cent      =',num2str(y(4)) ] );
disp( ['50.0 per cent      =',num2str(y(5)) ] );
disp( ['90.0 per cent      =',num2str(y(6)) ] );
disp( ['95.0 per cent      =',num2str(y(7)) ] );
disp( ['97.5 per cent      =',num2str(y(8)) ] );
disp( ['99.0 per cent      =',num2str(y(9)) ] );

% p   = [0.025 0.975 ];                      
% tmp = norminv(p,0,1)./sqrt(t);                       
% disp( ['95% confidence interval',num2str(tmp) ] );


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
% Panel (a)
subplot(2,2,1)
hist(rca,21)
title('(a) I[0] vs I[0]');
ylabel('$f(\widehat{\rho})$');
xlabel('$\widehat{\rho}$');
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
box off;

     
% %--------------------------------------------------------%
% % Panel (b)
subplot(2,2,2)
hist(rcb,21)
title('(b) I[1] vs I[1]');
ylabel('$f(\widehat{\rho})$');
xlabel('$\widehat{\rho}$');
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
box off;
  
%--------------------------------------------------------%
% Panel (c)
subplot(2,2,3)
hist(rcc,21)
title('(c) I[1] vs I[2]');
ylabel('$f(\widehat{\rho})$');
xlabel('$\widehat{\rho}$');
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
box off;

%--------------------------------------------------------%
% Panel (d)
subplot(2,2,4)
hist(rcd,21)
title('(d) I[2] vs I[2]');
ylabel('$f(\widehat{\rho})$');
xlabel('$\widehat{\rho}$');
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
box off;


 
% Print the tex file to the relevant directory
laprint(1,'spurious','options','factory');



