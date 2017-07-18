%========================================================================
%
%      Foreign exchange market efficient
%
%========================================================================

clear all;
clc;

% Load data
load sf.dat
f = sf(:,1);
s = sf(:,2);


% Spread
y = trimr(s,3,0)-trimr(f,0,3);

%********************************************************************
%***
%***     Generate graph
%***
%********************************************************************
dates = 1979+(4/12):1/12:2012-(4/12);

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

subplot(1,2,1)
plot(dates,trimr(s,3,0),'-k',dates,trimr(f,0,3),'--k')
title('Spot and Forward Rates');
ylim([1.0 2.8]) 
xlim([1978 2012]);
box off


subplot(1,2,2)
plot(dates,y,'-k');
title('Spread');
ylim([-0.4 0.3]);
xlim([1978 2012]);
box off
      
% Print the tex file to the relevant directory
% laprint(1,'fxefficiency','options','factory');

t = length(y)-3;
x = [ ones(t,1) trimr(y,0,3) ];
y = trimr(y,3,0);

xxi     = inv(x'*x);
bhat    = x\y;
uhat    = y-x*bhat;

% Plot first 6 autocorrelations of uhat and y(t-3)*uhat
figure(2)
autocorr(uhat,6,[ ],2)

figure(3)
autocorr(x(:,2).*uhat,6,[ ],2)


[~,cols] = size( x );
k        = 1;        
xuhat    = zeros(t,cols);
    
    for j = 1:cols
        
        xuhat(:,k) = x(:,j).*uhat;
        k=k+1;
    end      
s2      = uhat'*uhat/t;
seols   = sqrt(diag(s2*xxi));
white   = xxi*(xuhat'*xuhat)*xxi;
sewhite = sqrt(diag(white));

% Initial estimate of optimal lag length
P = floor(4*(t/100)^(2/9)); 
disp(['Initial estimate of lag length = ' num2str(P) ]);
J0 = xuhat'*xuhat
J1 = 0;
for j = 1:P
	Gj = trimr(xuhat,j,0)'*trimr(xuhat,0,j);
	J0 = J0 + Gj + Gj';
	J1 = J1 + 2*j*Gj;
end
J1
% Revised estimate of optimal lag length
i  = [ 1; ones(cols-1,1) ];
v0 = i'*J0*i
v1 = i'*J1*i
P  = floor(1.1447*((v1/v0)^2*t)^(1/3));
disp(['Revised estimate of lag length = ' num2str(P) ]);

% Compute Newey-West estimate of variance
JT = xuhat'*xuhat; 
for j = 1:P
	Gj = trimr(xuhat,j,0)'*trimr(xuhat,0,j);
	JT = JT + (1-j/(P+1))*(Gj + Gj');
end 

varNW = xxi*JT*xxi;
seNW = sqrt(diag(varNW));

disp( '    Coeff     se OLS    se White  se Newey-West ');
disp( [ bhat seols sewhite seNW ] );

% Wald test that both beta1 and beta2 = 0
wd = bhat'*inv(varNW)*bhat;
disp(['Wald statistic          = ',num2str(wd) ]); 
disp(['p-value                 = ',num2str(1-cdf('chi2',wd,2)) ]);

