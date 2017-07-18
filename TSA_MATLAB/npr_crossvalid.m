%=========================================================================
%
%     Cross validation approach to compute optimal bandwidth.
%
%=========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',1) );

% Simulate the model      
t  = 500;
ut = 0.1*randn(t,1);                %*     N(0,0.1^2)                      
xt = sort(-2 + 4*rand(t,1),1);      %*     U(-2,2)                         
mx = 0.3*exp( -4*(xt + 1).^2 ) + 0.7*exp( -16*(xt - 1).^2 );  
yt  = mx + ut;

% Compute weight function for cross validation (exclude lowest and highest 5% from sample)
xts = sort(xt,1);
xl  = xts( round(0.05*t));
xu  = xts( round(0.95*t));
wt  = zeros(t,1);
for i = 1:t
    if (xt(i)-xl>0)
        if (xt(i)-xu<0)
            wt(i)=1;
        end
    end
end

% Perform a grid search for h over the grid [0.01, 1.0]   
hopt = 0.001:0.001:0.3;
n    = length(hopt);
ace  = zeros(n,1);

for i = 1:n
    
    h   = hopt(i);
    tmp = (repmat(xt,1,t) - repmat(xt',t,1))/h; 
    fx  = mean( diagrv( normpdf(tmp)'/h,zeros(t,1) ) )';                % Leave-one-out   
   	fyx = mean( diagrv( normpdf(tmp)'.*repmat(yt,1,t)/h,zeros(t,1)) )'; % Leave-one-out   
    mx  = fyx./fx;

    ace(i) = sum( ((yt - mx).^2 ).*wt );

end


% Repeat calculations without cross-validation  
nce  = zeros(n,1);

for i = 1:n

    h   = hopt(i);
    tmp = (repmat(xt,1,t) - repmat(xt',t,1))/h; 
    fx  = mean( normpdf(tmp)'/h  )';                    
   	fyx = mean( normpdf(tmp)'.*repmat(yt,1,t)/h )';  
    mx  = fyx./fx;

    nce(i) = sum( ((yt - mx).^2 ).*wt );

end

[value,index] = min(ace);
disp(['\nBandwidth (with cross-validation)             = ', num2str(index) ] );
disp(['\nObjective function (with cross-validation)    = ', num2str(value) ] );
[value,index] = min(nce);
disp(['\nBandwidth (without cross-validation)          = ', num2str(index) ] );
disp(['\nObjective function (without cross-validation) = ', num2str(value) ] );

 
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

plot(hopt,ace,'-k',hopt,nce,'--k','LineWidth',0.75);
ylabel('$\mathcal{S}(h)$');
xlabel('$h$');
box off;
axis tight

%laprint(1,'crossvalidation','options','factory');

