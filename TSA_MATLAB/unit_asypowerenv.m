%=========================================================================
%
%   Approximate asymptotic power envelope and power curves 
%
%=========================================================================
function unit_asypowerenv( )

    % Critical values to compute the power envelope with cv=0 -> size of 0.05           
    cv   = seqa(-30,1,30);  

    n    = length( cv );
    t    = 1000;                       
    reps = 50000;                  
 
    % Detrending parameters
    x    = [ ones(t,1) seqa(1,1,t)' ];      
    cbar =-13.5;                     

    % Allocate memory
    pcc    = zeros( reps,n ); 
    pc0    = zeros( reps,n ); 
    pcbarc = zeros( reps,n ); 
    pcbar0 = zeros( reps,1 );
    dfols  = zeros( reps,n ); 
    dfols0 = zeros( reps,1 );
    dfgls  = zeros( reps,n ); 
    dfgls0 = zeros( reps,1 );
   
    h = waitbar(0,'Progress ....');
    for k = 1:n

        waitbar(k/n)
        RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234) );
        c = cv(k);  
        %disp( cv(k) )
    
        for j=1:reps

            u  = randn(t,1);
            yc = recserar(u,u(1),1+c/t);
            y0 = recserar(u,u(1),1);

            pcc(j,k)    = Ptest(yc,x,c,0);
            pc0(j,k)    = Ptest(y0,x,c,0);
            pcbarc(j,k) = Ptest(yc,x,cbar,0);
            dfols(j,k)  = trimr(adftests(glsdetrend(yc,x,-t),0),1,0);         
            dfgls(j,k)  = trimr(adftests(glsdetrend(yc,x,cbar),0),1,0);  
            
            if k==1
                
                pcbar0(j) = Ptest(y0,x,cbar,0); 
                dfols0(j) = trimr(adftests(glsdetrend(y0,x,-t),0),1,0);
                dfgls0(j) = trimr(adftests(glsdetrend(y0,x,cbar),0),1,0);
            end


        end

    end
    close(h);
    
    % Compute rejection frequencies for alternative detrending methods
    rejpc    = [ mean(bsxfun(@lt,pcc,quantile(pc0,0.05)))  0.05 ];
    rejpcbar = [ mean(bsxfun(@lt,pcbarc,quantile(pcbar0,0.05)))  0.05 ];
    rejdfols = [ mean(bsxfun(@lt,dfols,quantile(dfols0,0.05)))  0.05 ];
	rejdfgls = [ mean(bsxfun(@lt,dfgls,quantile(dfgls0,0.05)))  0.05 ];
    
    %**********************************************************************
    %***
    %***     Generate graph
    %***
    %**********************************************************************

    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;
    
    dvec = [ cv 0 ];
    
    %--------------------------------------------------------%
    % dvec,rejpc,'-k',     ...
    plot( dvec,rejpcbar,'--k', ...
         dvec,rejdfols,'-.k', ...
         dvec,rejdfgls,':k');
    ylabel('Power')
    ylim( [0.0 1.0 ] );
    set(gca,'ytick',[0.0 0.2 0.4 0.6 0.8 1]);
    xlim( [ -30 0 ] );
    set(gca,'xtick',[-30 -20  -10 0]);
    box off

    % Print the tex file to the relevant directory
    %laprint(1,'asypowerenv','options','factory');

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
%-------------------------------------------------------------------------
%  P-test
%-------------------------------------------------------------------------
function pt = Ptest(y,x,cbar,k)

    n  = length(y);
    uc = glsdetrend(y,x,cbar); 
    u0 = glsdetrend(y,x,0);
    s2 = ar1rvar(uc,k);
    uc = [ uc(1) ; trimr(uc,1,0)-(1+cbar/n)*trimr(uc,0,1) ];
    u0 = [ u0(1) ; trimr(u0,1,0)-trimr(u0,0,1) ];

    pt = (uc'*uc-(1+cbar/n)*u0'*u0)/s2;

end

%-------------------------------------------------------------------------
%  ar1rvar
%-------------------------------------------------------------------------
function s2 = ar1rvar(u,k)

    du = trimr(u,1,0)-trimr(u,0,1); 
    x  = trimr(u,0,1);

    if k>0 
        
        tmp = [x lagmatrix(du,1:1:k)];
        
        x = trimr( tmp,k,0 ); 
    end

    b = x\trimr(du,k,0); 
    e = trimr(du,k,0) - x*b; 
    s2 = e'*e/length(e);

    if k>0
        
        s2 = s2/(1-sum(trimr(b,1,0)))^2; 
    
    end    

end

