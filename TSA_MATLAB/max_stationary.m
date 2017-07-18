%=========================================================================
%
%   Maximum likelihood estimation of the stationary distribution of the 
%   CIR model of interest rates using Ait Sahalia's (1996) data.
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
    load eurodollar.mat
     % Load data: time series object 'tbr'
    %load FREDdata


    rt = sort( eurodata(:,4)*100 );
    %rt = tbr.Data;
    t  = length( rt );
    pstart = [ 5.6556    0.6763];

    options = optimset('LargeScale', 'off', 'Display', 'iter');     
    [phat,~,~,~,~,hessian] = fminunc(@(p) neglog(p,rt),pstart,options);
    %[phat] = fminsearch(@(p) neglog(p,rt),pstart,options);
    
    hessinv = inv(hessian);
    disp( 'Parameter estimates')
    disp( phat )
    disp( 'Standard errors based on inverse Hessian')
    disp( [sqrt( hessinv(1,1) ) sqrt( hessinv(2,2) )] )


    %********************************************************************
    %***
    %***     Generate graph
    %***
    %********************************************************************

    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;

    [fcdf,x] = ecdf(rt); 
    [f,bins] = ecdfhist(fcdf,x,51); 
    bar(bins,f,'hist');
    %[n,xout]=hist(rt,51);
    %bar(xout,n/t)
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k');
    box off
    hold on
    plot(rt,gampdf(rt,phat(1),1/phat(2)),'-k','LineWidth',0.75)
    ylabel('$f(r)$');
    xlabel('$r$');
    set(gca,'YTick',[] );
    axis tight
    hold off
    
    % Print the tex file to the relevant directory
    %laprint(1,'cirstat','options','factory');

    
end
%
%--------------------------- Functions -----------------------------------
% 

%-------------------------------------------------------------------------
% Likelihood wrapper function
%-------------------------------------------------------------------------
function f = neglog( p,data )

    f = -sum( lnlt( p,data) );
end

%-------------------------------------------------------------------------
% Likelihood function for stationary distribution of CIR model
%-------------------------------------------------------------------------
function f = lnlt(p,data)

    v = abs(p(1));
    w = abs(p(2));
    %f = log( gampdf(data,v,w ) );
    f = v*log( w ) - gammaln( v ) + (v-1)*log(data) - w*data;
end
