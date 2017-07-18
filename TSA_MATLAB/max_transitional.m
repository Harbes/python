%=========================================================================
%
%   Maximum likelihood estimation of the transitional distribution of the 
%   CIR model of interest rates using Ait Sahalia's (1996) data.
%
%=========================================================================

function max_transitional(  )

    clc
    clear all
    
    % Load data (5505x4 array called eurodata, 1 Jun 1973 - 25 Feb 1995)
    %   1. year
    %   2. day
    %   3. date stamp
    %   4. interest rates
    load eurodollar.mat
 
    dt = 1/250;
    rt = eurodata(:,4);

    % Starting values
    x  = rt(1:end-1);
    dx = diff(rt);           
    dx = dx./x.^0.5;
    regressors = [dt./x.^0.5, dt*x.^0.5];
    drift = regressors\dx;
    res = regressors*drift - dx;
    alpha = -drift(2);
    mu = -drift(1)/drift(2);
    sigma = sqrt(var(res, 1)/dt);
    
    p0 = [abs(alpha) abs(mu) abs(sigma)]; 
    
    % Estimation based on scaled Bessel function
    options = optimset('LargeScale', 'off', 'Display', 'iter');     
    [phat1,~,~,~,~,hessian] =  fminunc(@(p) cir1(p,rt),p0,options);  
    
    hessinv = inv(hessian);
    disp( 'Parameter estimates')
    disp( phat1 )
    disp( 'Standard errors based on inverse Hessian')
    disp( [sqrt( hessinv(1,1) ) sqrt( hessinv(2,2) ) sqrt( hessinv(3,3) )] )
 
    % Estimation based on ncx2pdf function
    
    options = optimset('LargeScale', 'off', 'Display', 'iter');     
    [phat,~,~,~,~,hessian] =  fminunc(@(p) cir2(p,rt),phat1,options);  
    
    hessinv = inv(hessian);
    disp( 'Parameter estimates')
    disp( phat )
    disp( 'Standard errors based on inverse Hessian')
    disp( [sqrt( hessinv(1,1) ) sqrt( hessinv(2,2) ) sqrt( hessinv(3,3) )] )
 
    
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
    tmp1 = (phat(3).^2)*tmp0; 
    
    plot(tmp0,tmp1,'-k','LineWidth',0.75);
    hold on
    scatter(rlag,rnow.^2,'.k')
    hold off
    ylabel('$r^2_{t+\Delta t}$');
    xlabel('$r_t$');
    set(gca,'YTick',[] );
    axis tight
    hold off
    
    % Print the tex file to the relevant directory
    laprint(1,'diffusion','options','factory');

end
    


%
%--------------------------- Functions -----------------------------------
% 


%-------------------------------------------------------------------------
% Likelihood function for transitional distribution of CIR model (Bessel)
%-------------------------------------------------------------------------
function f = cir1( p,data )

    alpha = abs(p(1));
    mu    = abs(p(2));
    sigma = abs(p(3));
    
    rnow = data(2:end);
    rlag = data(1:end-1);
    dt   = 1/250;
    
    c = 2*alpha/(sigma^2*(1-exp(-alpha*dt)));
    q = 2*alpha*mu/sigma^2-1;
    u = c*exp(-alpha*dt)*rlag;
    v = c*rnow;    
    
    lf = log(c)-u-v+0.5*q*log(v./u)+log(besseli(q,2*sqrt(u.*v),1))+2*sqrt(u.*v);
    f  = -sum( lf );
end

    
%-------------------------------------------------------------------------
% Likelihood function for transitional distribution of CIR model
% (Chi-square)
%-------------------------------------------------------------------------
function f = cir2( p,data )

    alpha = abs(p(1));
    mu    = abs(p(2));
    sigma = abs(p(3));
    
    rnow = data(2:end);
    rlag = data(1:end-1);
    dt   = 1/250;
 
    c  = 2*alpha/(sigma^2*(1-exp(-alpha*dt)));
    q  = 2*alpha*mu/sigma^2-1;
    u  = c*exp(-alpha*dt)*rlag;
    v  = c*rnow; 
    nc = 2*u;
    df = 2*q+2; 
    s  = 2*v;
    
    gpdf = ncx2pdf( s,df,nc );
    ppdf = 2*c*gpdf;

    f = -sum(log( ppdf ) );
end



