%=========================================================================
%
%   Maximum likelihood estimation of the transitional distribution of the 
%   Vasciek model of interest rates using Ait Sahalia's (1996) data.
%
%=========================================================================

function max_stationary( )
clear all;
clc;

% Load data (5505x4 array called eurodata, 1 Jun 1973 - 25 Feb 1995)
%   1. year
%   2. day
%   3. date stamp
%   4. interest rates
load eurodata.mat


rt = eurodata(:,4)*100;
Params = [ 0.436003 0.061221 0.149221 ];
f = -sum( lnlt(Params,rt) )

end

function f = lnlt(p,data)



    alpha = p(1);
    mu    = p(2);
    sigma = p(3);

    w= 2*alpha/sigma^2;
    v = 2*alpha*mu/sigma^2;
  
    f = -v*log( w ) - gammaln( v ) + (v-1)*data - w*data;
end
