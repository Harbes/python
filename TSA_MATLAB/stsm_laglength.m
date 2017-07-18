%=========================================================================
%
%   Simulate an AR(3) model and compute the optimal lag length
%
%==========================================================================



clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',1) );

% Parameters
NRep = 2000;
t    = 300;                       
pmax = 7;      
mu   = 0.0;
phi1 = 0.2;    
phi2 = -0.15;    
phi3 = 0.05;

acount = zeros( pmax,1 );
scount = zeros( pmax,1 );
hcount = zeros( pmax,1 );


for k = 1:NRep

    % Generate data
    yt = zeros( t,1 );

    for i = 4:t
    
        yt(i) = mu+phi1*yt(i-1)+phi2*yt(i-2)+phi3*yt(i-3)+sqrt(0.5)*randn;
    end
    
    % Set up lags and make sure that T is constant
    ylag = lagmatrix(yt,0:pmax);
    ylag = ylag(101:end,:);
    

    % Loop over the lags (yt is first column of ylag)
    aic = zeros( pmax,1 );
    sic = zeros( pmax,1 );
    hic = zeros( pmax,1 );
    y   = ylag(:,1);
    tt  = length( y );

    for j = 2:pmax+1

        x = [ ones(tt,1) ylag(:,2:j)];
        k = size( x,2 );
        b = x\y;
        e = y-x*b;
    
        aic(j-1) = log(e'*e/tt) + 2*k/tt;
        sic(j-1) = log(e'*e/tt) + k*log(tt)/tt;
        hic(j-1) = log(e'*e/tt) + 2*k*log(log(tt))/tt;
    end


    [~,ind]=min(aic);
    acount(ind) = acount(ind)+1;
    [~,ind]=min(sic);
    scount(ind) = scount(ind)+1;
    [~,ind]=min(hic);
    hcount(ind) = hcount(ind)+1;

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

bar([acount hcount scount])
colormap gray

xlabel('Lag Order');

laprint(1,'infocrit','options','factory');





