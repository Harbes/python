%=========================================================================
%
%   Approximate asymptotic power envelope and power curves of ADF test
%   using alternative detrending methods
%
%
%=========================================================================
function unit_asypower( )

    clear all;
    clc;

    % Critical values to compute the power envelope with cv=0 -> size of 0.05           
    cv   = seqa(-30,1,31);  

    n    = length( cv );
    t    = 1000;                       
    reps = 50000;                  
 
    % Detrending parameters
    x    = [ ones(t,1) seqa(1,1,t)' ];      
    cbar =-13.5;                     

    dfols  = zeros( reps,n ); 
    dfgls  = zeros( reps,n ); 

    h  = waitbar(0,'Progress...');
    
    for k = 1:n

        
        RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234) );
        c = cv(k);  
        
    
        for j=1:reps

            u  = randn(t,1);
            yc = recserar(u,u(1),1+c/t);

            dfols(j,k) = trimr(adftests(glsdetrend(yc,x,-t),0),1,0);         
            dfgls(j,k) = trimr(adftests(glsdetrend(yc,x,cbar),0),1,0);       

        end
        waitbar(k / n);

    end
    close(h)
    % Compute rejection frequencies for alternative detrending methods
    rejdfols = mean(bsxfun(@lt,dfols,quantile(dfols(:,end),0.05)));
	rejdfgls = mean(bsxfun(@lt,dfgls,quantile(dfgls(:,end),0.05)));

    
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
    plot(cv,rejdfols,'-k', ...
         cv,rejdfgls,'--k');
    ylabel('Power')
    ylim( [0.0 1.0 ] );
    set(gca,'ytick',[0.0 0.2 0.4 0.6 0.8 1]);
    xlim( [ -30 0 ] );
    set(gca,'xtick',[-30 -20  -10 0]);

    box off

    % Print the tex file to the relevant directory
    %laprint(1,'asypower','options','factory');
end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Detrending function: 
%       cbar = -7 constant; 
%       cbar = -13.5 linear trend  
%       cbar = -T for OLS detrending
%-------------------------------------------------------------------------
function u = glsdetrend( y,x,cbar )

    t = length( y );
    yc = [ y(1)  ; (trimr(y,1,0)-(1+cbar/t)*trimr(y,0,1)) ];
    xc = [x(1,:) ; (trimr(x,1,0)-(1+cbar/t)*trimr(x,0,1)) ];
    
    b  = xc\yc;
    u  = y - x*b;
end
%-------------------------------------------------------------------------
%  ADF coefficient and t tests: u must be already de-trended
%-------------------------------------------------------------------------
function adfstats= adftests(u,k)


    du = trimr(u,1,0) - trimr(u,0,1);
    x  = trimr(u,0,1);

    % Set up lags of du
    if k>0;
    
        ldu = lagmatrix(du,1:1:k);
        x   = trimr( [x ldu],k,0 );
    end
    
    xxi = inv(x'*x); 
    b   = xxi*x'*trimr(du,k,0); 
    e   = trimr(du,k,0)-x*b; 
    s2  = e'*e/length(e);

    adfstats = [length(u)*b(1) ; b(1)/sqrt(s2*xxi(1,1))];

end
