%==========================================================================
%
%     Vuoung's Nonnested Test
%
%==========================================================================
function nls_money(  )

    clear all;
    clc;

    % Load data
    load moneydemand.mat
    
    mt = m2./cpi;
    yt = gdp./cpi;
    rt = tbill/100;
    t  = length( mt );       	% Define the sample size      

    % Estimate model 1
    y    = mt;                  	% Estimate Model 1    
    x    = [ones(t,1) rt yt];
    b    = x\y;
    sig2 = ((y - x*b)'*(y - x*b))/t;
    lf1  = model1( [ b; sig2],y,x );

    % Estimate Model 2 
    y    = log(mt);              	   
    x    = [ones(t,1) log(rt) log(yt)];
    b    = x\y;
    sig2 = ((y - x*b)'*(y - x*b))/t;
    lf2  = model2( [ b; sig2],y,x );

    % Perform Vuong's test
    dt   = lf1 - lf2;
    dbar = mean( dt );
    sbar = std( dt,1 );

    V = sqrt(t)*(dbar/sbar);
    
    disp('Vuongs test and p-value');
    disp('------------------------');
    disp([ V  normcdf(V,0,1) ] );
    
end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%   Likelihood function model 1
%-------------------------------------------------------------------------

function lf = model1( b,y,x )

    u  = y-b(1)*x(:,1)-b(2)*x(:,2)-b(3)*x(:,3);
    lf = -0.5*log(2*pi) - 0.5*log( b(4) ) - u.^2/(2*b(4)) + log(abs(1.0));
    
end
%-------------------------------------------------------------------------
%   Likelihood function model 2
%-------------------------------------------------------------------------
function lf = model2( b,y,x )

    u  = y - b(1)*x(:,1) - b(2)*x(:,2) - b(3)*x(:,3);
    lf = -0.5*log(2*pi) - 0.5*log( b(4) ) - u.^2/(2*b(4)) ...
        + log( abs( 1.0./exp( y ) ) );
    
end




