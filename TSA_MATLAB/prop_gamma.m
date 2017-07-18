%=========================================================================
%
%    Program to demonstrate the Lindberg-Feller central limit theorem
%    using a regression model with gamma distributed errors.
%
%=========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',123456) )

t      = 500;
NDraws = 5000;                       % Number of draws                                    


b0    = [1; 2];                      % Population parameters   
rho   = 0.25;                        % Parameters of the gamma distribution     
alpha = 0.1;


trigs = randg( rho,[t NDraws]);
x0    = [ones(t,1) randn(t,1)];

% For the sample size of T = 100        
x  = x0(1:100,:);
z1 = zeros(NDraws,1);
z2 = zeros(NDraws,1);

for i = 1:NDraws;

    u    = alpha*trigs(1:100,i) - rho*alpha;  
    y    = x*b0 + u;
    b    = x\y;
    e    = y - x*b;
    s2   = e'*e/100;
    vcov = s2*inv(x'*x);
    z1(i)= (b(1) - b0(1))/sqrt(vcov(1,1));
    z2(i)= (b(2) - b0(2))/sqrt(vcov(2,2));

end
     
% For the sample size of T = 500        
x  = x0(1:500,:);
z3 = zeros(NDraws,1);
z4 = zeros(NDraws,1);

for i = 1:NDraws;

    u    = alpha*trigs(1:500,i) - rho*alpha;  
    y    = x*b0 + u;
    b    = x\y;
    e    = y - x*b;
    s2   = e'*e/500;
    vcov = s2*inv(x'*x);
    z3(i)= (b(1) - b0(1))/sqrt(vcov(1,1));
    z4(i)= (b(2) - b0(2))/sqrt(vcov(2,2));

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

subplot(2,2,1)
hist(z1,21);
title('(a) Distribution of $z_{\widehat{\beta}_0}$ (T = 100)');
ylabel('$f(z_{\widehat{\beta}_0})$');
%xlabel('$Z_{\widehat{\beta}_0}$');
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')

subplot(2,2,2)
hist(z2,21);
title('(b) Distribution of $z_{\widehat{\beta}_1}$ (T = 100)');
ylabel('$f(z_{\widehat{\beta}_1})$');
%xlabel('$Z_{\widehat{\beta}_1}$');
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')

subplot(2,2,3)
hist(z3,21);
title('(c) Distribution of $z_{\widehat{\beta}_0}$ (T = 500)');
ylabel('$f(z_{\widehat{\beta}_0})$');
%xlabel('$Z_{\widehat{\beta}_0}$');
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')

subplot(2,2,4)
hist(z4,21);
title('(d) Distribution of $z_{\widehat{\beta}_0}$ (T = 500)');
ylabel('$f(z_{\widehat{\beta}_0})$');
%xlabel('$Z_{\widehat{\beta}_0}$');
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')


% Print the TEX file
laprint(1,'cltreg','options','factory');


