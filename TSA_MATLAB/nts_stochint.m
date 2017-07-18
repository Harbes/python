%=========================================================================
%
%   Program to simulate the distribution of a stochastic integral
%
%=========================================================================

clear all;
clc;

RandStream.setGlobalStream( RandStream('mt19937ar','seed',15) )

imoment = 1;

t     = 500;
nreps = 10000;      
delta = 0.0;
phi   = 1.0;
sig2  = 1.0;
sig   = sqrt(sig2);
y0    = 0.0;


m = zeros(nreps,1);        

for i = 1:nreps

    % Generate random walk and discard first element
    v  = sqrt(sig2)*randn(t+1,1);
    y  = trimr( recserar(delta + v , y0 , phi),1,0);    
    v  = trimr(v,1,0);                                      
    if imoment == 1;                                                                                 

        % Sample mean with standardization based on t^(-0.5)
        m(i) = (1/sig^2)*sum( trimr(y,0,1).*trimr(v,1,0))*t^(-1);                                  

    elseif imoment == 2; 
        
        % Standardized sum of y(t-1)^2 * v(t) with standardization based on sig^(-3)xt^(-3/2)                         
        m(i) = (1/sig^3)*sum( trimr(y.^2,0,1).*trimr(v,1,0))*t^(-3/2);     
    end
end

disp( ['Sample size = ',num2str(t) ] );
disp(' ' );



if imoment == 1;

    xi = seqa(-0.99,0.1,201);

    
    fhat  = ksdensity( m,xi );
    %ftrue = chi2pdf(xi,1);
    tmp = 2*xi + 1;                                                    
    jacobian = 2;
    ftrue = jacobian*(tmp.^(-0.5)).*exp(-tmp/2)/(gamma(0.5)*sqrt(2));


    disp( ['Sample mean of m           = ',num2str(mean(m)) ] );
    disp( ['Theoretical mean of m      = ',num2str(0.0) ] );
    disp( ['Sample variance of m       = ',num2str(std(m)^2) ] );
    disp( ['Theoretical variance of m  = ',num2str(0.5) ] );


elseif imoment == 2; 

    xi = seqa(-5,0.1,101);

    fhat = ksdensity( m,xi );
    ftrue = normpdfn(xi);

    disp( ['Sample mean of m           = ',num2str(mean(m)) ] );
    disp( ['Sample variance of m       = ',num2str(std(m)^2) ] );



end

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

plot(xi,ftrue,'-k', ...
     xi,fhat,'--k');
box off
ylim([0 2]);
xlim( [-1 4]);
set(gca,'xtick',[-1 0 1 2 3 4]);
set(gca,'ytick',[0.0 0.5 1.0 1.5 2.0]);
xlabel('$m$');
ylabel('$f(m)$');

% Print the tex file to the relevant directory
%laprint(1,'stochint','options','factory');


