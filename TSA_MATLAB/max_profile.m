%=========================================================================
%
%     Program to plot the profile likelihood for the portfolio
%     diversification model
%     
%
%=========================================================================


clear all;
clc;

% Asset price data from 6 August 2010 to 2 January 2001 
%(note that the data are in reverse order ie from recent to past	

%load data

load diversify.mat;

% Select appropriate sample  
pt_apple     = pt_apple(1:2413);
pt_ford      = pt_ford(1:2413);

%     Compute percentage returns     
r_apple = 100*( log(pt_apple(2:end,:)) - log(pt_apple(1:(end-1),:)) );

r_ford  = 100*( log(pt_ford(2:end,:)) - log(pt_ford(1:(end-1),:)) );


%     Compute maximum likelihood estimates    

m1 = mean(r_apple);

m2 = mean(r_ford);

s11 = mean((r_apple - mean(r_apple)).^2);

s22 = mean((r_ford - mean(r_ford)).^2);

c   = corrcoef(r_apple,r_ford);

rho   = c(1,2);
disp(['Value of rho = ' num2str(rho)]);

%     Generate the profile log-likelihood     

y = [r_apple, r_ford];

r = ( -0.9 :0.01:-0.9+0.01*(186-1) )';

a = zeros(length(r),1);

mean = [ m1, m2];

for i = 1:length(r)
    
    s12 = r(i)*sqrt(s11)*sqrt(s22);      % covariance between apple and ford
    covariance = [s11, s12; s12, s22];   % s12 = s21
    a(i) = -ecmnobj(y, mean, covariance)/length(y);  % Compute average log-likelihood value for alternative values of rho

end

%     Plot the profile log-likelihood        

%********************************************************************
%***
%***     Generate graph
%***
%********************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

plot(r,a,'k');

xlabel('theta1');
ylabel('Average lnl');

