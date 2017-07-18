%=========================================================================
%
%   Program to demonstrate the asymptotic normality of the QMLE 
%   where the true distribution is Student t
%
%=========================================================================

clear all;
clc;
format short;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234) )

mu  = 10;      
sig = 1;       
gam = 5;       

theta0 = [mu ; sig^2 ; gam;  ];

NSim = 50000;
NObs = [ 50, 100, 200, 400, 800, 1600 ];

% Allocate vectors to hold standardized values of mu and sig^2
z1 = zeros( NSim,length(NObs) );
z2 = zeros( NSim,length(NObs) );          
    
for k = 1:1:length( NObs );
    
    t      = NObs(k);
    theta1 = zeros( NSim,2);

    for n = 1:1:NSim;

        % Generate data from t distribution
        v = trnd( gam,t,1 );                                
        y = mu + sig*sqrt((gam-2)/gam)*v;                  

        % Estimate parameters of misspecified model
        m  = mean(y);                                   
        s2 = mean( (y - m).^2 );

        % Gradients of misspecified model
        g1 = (y - m)/s2;        
        g2 = -0.5/s2 + 0.5*(y - m).^2/s2^2;
        G  = [g1  g2];

        % OPG
        J = G'*G;
        
        % Hessian
        H = zeros(2,2);         
        
        H(1,1) = -t/s2;
        H(1,2) = -sum(y - m)/s2^2;
        H(2,1) = H(1,2);
        H(2,2) = t*0.5/s2^2 - sum((y - m).^2)/s2^3;

        % Information matrix
        I = -H;

        % QMLE standard errors
        se = sqrt( diag( inv(I)*J*inv(I) ) );
           
        theta1(n,:) = [m  s2];

        z1(:,k) = (theta1(:,1) - theta0(1))./se(1);     
        z2(:,k) = (theta1(:,2) - theta0(2))./se(2);    

    end;
            
    tmp = bsxfun(@minus,theta1,theta0(1:2)');
    rmse = sqrt(mean(tmp).^2);

    disp( ['Sample size of T    =' num2str(t)]);
    disp( ['True                =' num2str(theta0(1:2)')]);
    disp( ['Mean                =' num2str(mean(theta1))]);
    disp( ['Standard deviation  =' num2str(std(theta1))]);
    disp( ['RMSE                =' num2str(rmse)]);
    disp( ' ' );
    disp( '----------------------------------------');
    disp( ' ' );

end;


% Estimate the nonparametric density

zgrid = -5:0.1:5;
f1    = zeros( length(zgrid),length( NObs ) );
f2    = zeros( length(zgrid),length( NObs ) );

for k = 1:1:length( NObs )
    
    f1(:,k) = ksdensity( z1(:,k),zgrid );
    f2(:,k) = ksdensity( z2(:,k),zgrid );
    
end

%**********************************************************************
%***
%***     Generate graphs
%***
%**********************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;


%--------------------------------------------------------%
% Panel (a)
subplot(1,2,1);
plot(zgrid,normpdf(zgrid),'-k', ...
     zgrid,f1(:,1),'-.k', ...
     zgrid,f1(:,3),'--k');
title('(a) $f(z_1)$ ');
ylabel('$f(z_1)$');
xlabel('$z_1$');
box off;

%--------------------------------------------------------%
% Panel (b)
subplot(1,2,2);
plot(zgrid,normpdf(zgrid),'-k', ...
     zgrid,f2(:,1),'-.k', ...
     zgrid,f2(:,3),'--k');
title('(b) $f(z_2)$ ');
ylabel('$f(z_2)$');
xlabel('$z_2$');
box off;


laprint(1,'qmleclt','options','factory');


