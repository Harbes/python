%=========================================================================
%
%   Program to compute the maximum likelihood estimates of the portfolio
%   diversification model
%
%   Asset price data from 6 August 2010 to 2 January 2001 (note that the
%   data are in reverse order ie from recent to past)
%
%=========================================================================

clear all
clc

% Load data

load diversify.mat

% Select appropriate sample  
pt_apple     = pt_apple(1:2413);
pt_ford      = pt_ford(1:2413);


% Compute percentage returns  
r_apple = 100*diff(log(pt_apple)); 
r_ford  = 100*diff(log(pt_ford));


% Compute statistics  
m1 = mean(r_apple);
m2 = mean(r_ford);

s11 = mean((r_apple - mean(r_apple)).^2);
s22 = mean((r_ford - mean(r_ford)).^2);
c   = corrcoef(r_apple,r_ford);
r   = c(1,2);


disp(['Sample mean (Apple)         = ' num2str(m1) ]);
disp(['Sample mean (Ford)          = ' num2str(m2) ]);

disp(['Sample variance (Apple)     = ' num2str(s11) ]);
disp(['Sample variance (Ford)      = ' num2str(s22) ]);

disp(['Sample correlation          = ' num2str(r) ]);


% Compute weights and risk of optimal portfolio

w1      = ( s22 - r*sqrt(s11*s22) ) / ( s11 + s22 - 2*r*sqrt(s11*s22) );
w2      = 1 - w1;
s2_port = w1^2*s11 + w2^2*s22 + 2*w1*w2*r*sqrt(s11*s22);

disp(['Optimal weight (Apple)      = ' num2str(w1) ]);
disp(['Optimal weight (Ford)       = ' num2str(w2) ]);

disp(['Risk of optimal portfolio   = ' num2str(s2_port) ]);


%********************************************************************
%***
%***     Generate graph
%***
%********************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

scatter(r_ford,r_apple,'.','k')

xlabel('Ford');
ylabel('Apple');

% Print TEX file
laprint(1,'diversify','options','factory');

