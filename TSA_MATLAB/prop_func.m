%=========================================================================
%
%     Program to demonstrate the convergence of the sample log-likelihood
%     function
%
%=========================================================================
clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',123456) )


t    = [5 20 100];
mu   = -3.0:0.1:3;
true = -0.5*(1 + mu.^2);
samp = zeros(length(mu),length(t));

for i = 1:length(t)
    
    yt   = randn(t(i),1);
    
    for j = 1:length(mu)
        samp(j,i) = -0.5*mean((yt-mu(j)).^2);
    end
end

%**************************************************************************
%**     
%**     Generate graph       
%**                                         
%**************************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
f1=figure(1);
clf;

plot(mu,true,'-k');
hold on; 
plot(mu,samp(:,1),'--k');
plot(mu,samp(:,2),'-.k');
plot(mu,samp(:,3),':k','LineWidth',0.75);
hold off;
ylabel('$\ln L(\mu)$');
xlabel('$\mu$');
legend('True','t=5','t=20','t=100','Location','South');
legend boxoff;
box off;
legend2latex(f1)
%axis tight;


%laprint(1,'norm','options','factory');



