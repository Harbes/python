%=========================================================================
%
%   Program to demonstrate the Law of Large Numbers 
%   (Exponential distribution example) 
%
%=========================================================================

clear all
clc

RandStream.setDefaultStream( RandStream('mt19937ar','seed',12) )

mu = 5;      
t  = 500;

% Generate exponential random numbers from a Gamma distribution
e = randg(1,[t 1]);
y = mu.*e;           

% Generate sample means from sample of size t=1,2,3,...,tmax    
ybar = zeros(t,1);

for i = 1:t 

    ybar(i) = mean(y(1:i));            

end

%**************************************************************************
%**
%**     Generate graph
%**
%**************************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;



tt = (1:1:t)';
plot(tt,mu*ones(t,1),'-k',...
     tt,ybar,'-k',...
     tt,4.80*ones(t,1),':k',...
     tt,5.20*ones(t,1),':k',...
    'LineWidth',0.75);

ylabel('$\overline{y}_T$');
xlabel('T');
set(gca,'XTick',0:100:500);
axis( [0 500 3 7 ]);
set(gca,'YTick',3:1:7);

%laprint(1,'wln','options','factory');


