%=========================================================================
%
%   Effect of the initial condition on mztOLS and mztGLS 
%   Union of rejections included
%
%=========================================================================

function unit_power0UR( )

    clear all;
    clc;
    
    % Alternative initial conditions
    u0v = [ 0 2.5 5 ];

    % Critical values to compute the power envelope with cv=0 -> size of 0.05           
    cv   = seqa(-30,1,31);  

    n    = length( cv );
    t    = 200;                       
    reps = 100000;                  
 
    % Test parameters
    x    = [ ones(t,1) seqa(1,1,t)' ];      
    cbar =-13.5; 
    cvgls = -2.867;
    cvols = -3.130;
    tau   = 1.038;


    mztols = zeros( reps,n,3 ); 
    mztgls = zeros( reps,n,3 ); 
    rejols = zeros( 3,n );
    rejgls = zeros( 3,n );
    rejur  = zeros( 3,n );
    
    for i = 1:length(u0v);
        disp( i );
        
        for k = 1:n

            RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234) )
            c = cv(k);     
            
            for j=1:reps

                u = [ u0v(i) ; randn(t,1) ];
                y = trimr( recserar(u,u(1),1+c/t),1,0 );

                mztols(j,k,i) = trimr( mtests(glsdetrend(y,x,-t),0),2,0 );         
                mztgls(j,k,i) = trimr( mtests(glsdetrend(y,x,cbar),0),2,0 );       

            end

        end
                
        % Compute rejection frequencies 
        rejols(i,:) = mean(bsxfun(@lt,mztols(:,:,i),cvols));
        rejgls(i,:) = mean(bsxfun(@lt,mztgls(:,:,i),cvgls));
        rejur(i,:)  = mean(bsxfun(@lt,mztols(:,:,i),cvols*tau) ...
                      | bsxfun(@lt,mztgls(:,:,i),tau*cvgls));


    end


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
    subplot(1,3,1)
    plot(cv,rejols(1,:),'-k', ...
         cv,rejgls(1,:),'--k',...
         cv,rejur(1,:),'-.k');
    title('(a) $v_0=0.0$')
    ylabel('Power')
    xlabel('c')
    ylim( [0.0 1.0 ] );
    set(gca,'ytick',[0.0 0.2 0.4 0.6 0.8 1]);
    xlim( [ -30 0 ] );
    set(gca,'xtick',[-30 -20  -10 0]);
    box off

    %--------------------------------------------------------%
    subplot(1,3,2)
    plot(cv,rejols(2,:),'-k', ...
         cv,rejgls(2,:),'--k',...
         cv,rejur(2,:),'-.k');
    title('(b) $v_0=2.5$')
     ylabel('Power')
    xlabel('c')
    ylim( [0.0 1.0 ] );
    set(gca,'ytick',[0.0 0.2 0.4 0.6 0.8 1]);
    xlim( [ -30 0 ] );
    set(gca,'xtick',[-30 -20  -10 0]);
    box off
    
    %--------------------------------------------------------%
    subplot(1,3,3)
    plot(cv,rejols(3,:),'-k', ...
         cv,rejgls(3,:),'--k',...
         cv,rejur(3,:),'-.k');
    title('(c) $v_0=5.0$ ')
     ylabel('Power')
    xlabel('c')
    ylim( [0.0 1.0 ] );
    set(gca,'ytick',[0.0 0.2 0.4 0.6 0.8 1]);
    xlim( [ -30 0 ] );
    set(gca,'xtick',[-30 -20  -10 0]);
    box off
    % Print the tex file to the relevant directory
    laprint(1,'unionrej','options','factory');

    
    
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

%-------------------------------------------------------------------------
%  M tests
%-------------------------------------------------------------------------
function mz = mtests( u,k )

    t  = length(u);
    s2 = ar1rvar( u,k );
    u2 = sum( u(1:t-1).^2/t^2 );

    mza = (u(t)^2/t-s2)/(2*u2);
    msb = sqrt(u2/s2);
    mzt = msb*mza;

    mz = [ mza ; msb; mzt ];

end

