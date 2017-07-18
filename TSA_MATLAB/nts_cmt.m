%==========================================================================
%
%   Properties of the Functional Central Limit Theorem for various moments
%
%==========================================================================
clear all
clc

RandStream.setDefaultStream( RandStream('mt19937ar','seed',15) );

ndraws  = 50000;
imoment = 1;             
t       = 500;

% Parameters
delta = 0.0;
phi   = 1.0;
sig2  = 1.0;
sig   = sqrt(sig2);
y0    = 0.0;

m = zeros(ndraws,1);         

for i = 1:ndraws

    % Random walk with first observation discarded
    y  = trimr( recserar(delta + sqrt(sig2)*randn(t+1,1),y0,phi) , 1 , 0);   

    if imoment == 1                                                                                

        % Sample mean with standardization based on t^(-0.5)
        m(i) = (1/sig)*sum(y)*t^(-3/2); 

    elseif imoment == 2

        m(i) = (1/sig)*sum( seqa(1,1,t)'.*y)*t^(-5/2); 

    end
end

disp(['Sample size =              = ', num2str(t) ])
disp(' ')
disp(['Sample mean of m           = ', num2str(mean(m)) ])
disp(['Theoretical mean of m      = ', num2str(0.0) ])
disp(['Sample variance of m       = ', num2str(std(m)^2) ])
disp(['Theoretical variance of m  = ', num2str(1/3) ])

x = seqa(-5,0.1,101);

%h = 1.06*std(m)/(t^0.2);
%h    = 1.06*minc(std(m)|((quantile(m,0.75)-quantile(m,0.25))/1.34))/(t^0.2);
%z    = bsxfun(@minus,m,x)/h;
%fhat = mean( 0.75*(1 - z^2).*(abs(z)<1) )/h;

fhat = ksdensity(m,x);
fnorm = normpdf(x);
plot(x,fnorm,x,fhat);


