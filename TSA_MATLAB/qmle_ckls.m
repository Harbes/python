%=========================================================================
%
%   Maximum likelihood estimation of the transitional distribution of the 
%   CKLS model of interest rates using Ait Sahalia's (1996) data.
%
%=========================================================================

function qmle_ckls(  )

    clc
    clear all
    
    % Load data (5505x4 array called eurodata, 1 Jun 1973 - 25 Feb 1995)
    %   1. year
    %   2. day
    %   3. date stamp
    %   4. interest rates
    load eurod.mat
 
    dt = 1/250;
    rt = eurodata(:,4);

    % Starting values
    x  = rt(1:end-1);
    dx = diff(rt);           
    dx = dx./x.^0.5;
    regressors = [dt./x.^0.5, dt*x.^0.5];
    drift = regressors\dx;
    res   = regressors*drift - dx;
    alpha = -drift(2);
    mu    = -drift(1)/drift(2);
    sigma = sqrt(var(res, 1)/dt);
    
    p0 = [abs(alpha) abs(mu) abs(sigma) 0.5]; 
    
    % Estimation based on scaled Bessel function
    options = optimset('LargeScale', 'off', 'Display', 'iter');     
    [phat,~,~,~,~,hessian] =  fminunc(@(p) ckls(p,rt),p0,options);  
    
    hessinv = inv(hessian);
    disp( 'Parameter estimates')
    disp( phat )
    disp( 'Standard errors based on inverse Hessian')
    disp( sqrt( diag(hessinv) ) )
 
 
    
    %********************************************************************
    %***
    %***     Generate graph
    %***
    %********************************************************************

    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;

    rlag = rt(1:end-1);
    rnow = rt(2:end);
    tmp0 = 0.03:0.01:0.25;
    tmp1 = phat(3)*tmp0.^phat(4);
   
    
    plot(tmp0,tmp1.^2,'-k','LineWidth',0.75);
    hold on
    scatter(rlag,(rnow-0.029).^2,'.k')
    hold off
    ylabel('$r^2_{t}$');
    xlabel('$r_{t-1}$');
    set(gca,'YTick',[] );
    axis tight
    hold off
%     
%     % Print the tex file to the relevant directory
%     laprint(1,'cklsdiffusion','options','factory');

end
    


%
%--------------------------- Functions -----------------------------------
% 

%-------------------------------------------------------------------------
% Quasi-likelihood function for transitional distribution of CKLS model  
%-------------------------------------------------------------------------
function lf = ckls( p,data )

    t = length(data)-1;
    
    alpha = abs(p(1));
    mu    = abs(p(2));
    sigma = abs(p(3));
    gamma = abs(p(4));
    
    rnow  = data(2:end);
    rlag  = data(1:end-1);
    dt    = 1/250;
    drift = alpha*(mu-rlag)*dt;
    stdev = sigma*rlag.^(gamma)*sqrt(dt);

    ut  =  (rnow - rlag - drift)./stdev;

    tmp = -0.5*ut.^2 - 0.5*log(stdev.^2) - 0.5*log(2*pi);
    lf  = -sum( tmp );
    
end

    
