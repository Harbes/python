%=========================================================================
%
%      Program to demonstrate a few nonlinear features  
%
%=========================================================================

function nlm_features( )

    clear all
    clc
    
    % Stationary time series
    %---------------------------------------------------------------------
    t = 40;
    y = 0.5 + zeros( t,1 );

    for i = 3:t

    y(i) = 0.9 + 0.7*y(i-1) - 0.6*y(i-2);

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

    subplot(1,2,1)
    plot(1:1:t,y,'-k','LineWidth',0.75);
    title( '(a)' )
    ylabel('$y_t$');
    xlabel('t');
    axis tight
    box off
	%xtics(0,40,10,0);
	%ytics(0,1.5,0.5,0);
    
     
    subplot(1,2,2)
    plot(trimr(y,0,1),trimr(y,1,0),'-k','LineWidth',0.75);
    title('(b)');
	xlabel('$y_{t-1}$');
	ylabel('$y_t$');
	%xtics(0,1.5,0.5,0);
	%ytics(0,1.5,0.5,0);
    box off
    

    
    % Print the tex file to the relevant directory
    %laprint(1,'nonlinar2','options','factory');


    % Limit Cycles (taken from Tong (1983, p.85)  
    % ---------------------------------------------------------------------
    t = 100;
    y = 2 + zeros(t,1);

    for i = 7:t

        if y(i-2) <= 3.05

            y(i) = 0.8023 + 1.0676*y(i-1) - 0.2099*y(i-2) + 0.1712*y(i-3) ...
                 - 0.4528*y(i-4) + 0.2237*y(i-5) - 0.0331*y(i-6);

        else

        y(i) = 2.2964 + 1.4246*y(i-1) - 1.0795*y(i-2) - 0.090*y(i-3);
        
        end
        

    end

    %**********************************************************************
    %***
    %***     Generate graph
    %***
    %**********************************************************************

    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(2);
    clf;

    subplot(1,2,1)
    plot(1:1:t,y,'-k','LineWidth',0.75);
    title( '(a)' )
    ylabel('$y_t$');
    xlabel('t');
    axis tight
    box off
	%xtics(0,100,20,0);
	%ytics(2,4,0.5,0);
    
     
    subplot(1,2,2)
    plot(trimr(y,0,1),trimr(y,1,0),'-k','LineWidth',0.75);
    title('(b)');
	xlabel('$y_{t-1}$');
	ylabel('$y_t$');
	%xtics(0,1.5,0.5,0);
	%ytics(0,1.5,0.5,0);
    box off
    
   
    % Print the tex file to the relevant directory
    %laprint(2,'limit','options','factory');

    
    % Strange attactor and chaos based on the Kaldor model 
    % (Lorenz (1989, p.130)) 
    %--------------------------------------------------------------------
    nobs = 10000;

    a = 5.0;
    c = 20.0;
    d = 0.01;
    e = 0.05;
    f = 280.0;
    g = 4.5;
    s = 0.21;
    
    alfa = 20.0;
    delta = 0.05;
    eta   = 0.00001;
    

    t = nobs + 5000;
    y = zeros(nobs,1);
    k = zeros(nobs,1);

    y(1) = 65.0;
    k(1) = 265.0;

    for i=1:t-1;
        invest = c*2^(-1/((d*y(i) + eta)^2)) + e*y(i) + a*(f/k(i))^g;
        saving = s*y(i);
        y(i+1) = y(i) + alfa*(invest - saving);
        k(i+1) = k(i) + invest - delta*k(i);
    end
    y = trimr(y,5000,0);
    k = trimr(k,5000,0);
    
    
    %**********************************************************************
    %***
    %***     Generate graph
    %***
    %**********************************************************************

    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(3);
    clf;

    subplot(2,2,1)
    plot(y,k,'.k','Markersize',1);
    title( '(a)' )
    ylabel('$k_t$');
    xlabel('$y_t$');
    axis tight
    box off
	%xtics(0,40,10,0);
	%ytics(0,1.5,0.5,0);
    
     
    subplot(2,2,2)
    plot(trimr(y,0,1),trimr(y,1,0),'.k','Markersize',1);
    title('(b)');
	xlabel('$y_{t}$');
	ylabel('$y_{t+1}$');
	%xtics(0,1.5,0.5,0);
	%ytics(0,1.5,0.5,0);
    box off
    
    subplot(2,2,3)
    plot(trimr(k,0,1),trimr(k,1,0),'.k','Markersize',1);
    title('(c)');
	xlabel('$k_{t}$');
	ylabel('$k_{t+1}$');
	%xtics(0,1.5,0.5,0);
	%ytics(0,1.5,0.5,0);
    box off
  
    % Show dependence on initial conditions
    nobs = 100;
    
    y0 = zeros( nobs,1 );
    k0 = zeros( nobs,1 );
    y0(1) = 65.0;
    k0(1) = 265.0;

    for i=1:nobs-1;
        invest = c*2^(-1/((d*y0(i) + eta)^2)) + e*y0(i) + a*(f/k0(i))^g;
        saving = s*y0(i);
        y0(i+1) = y0(i) + alfa*(invest - saving);
        k0(i+1) = k0(i) + invest - delta*k0(i);
    end
  
    y1 = zeros( nobs,1 );
    k1 = zeros( nobs,1 );
    y1(1) = 65.1;
    k1(1) = 265.1;

    for i=1:nobs-1;
        invest = c*2^(-1/((d*y1(i) + eta)^2)) + e*y1(i) + a*(f/k1(i))^g;
        saving = s*y1(i);
        y1(i+1) = y1(i) + alfa*(invest - saving);
        k1(i+1) = k1(i) + invest - delta*k1(i);
    end
  
    subplot(2,2,4)
    plot(1:1:nobs,y0,'-k',1:1:nobs,y1,'--k','LineWidth',0.75);
    title('(d)');
	xlabel('Time');
	ylabel('$y_{t}$');
    box off

    % Print the tex file to the relevant directory
    %laprint(3,'kaldor','options','factory');

    
    % Multiple equilibria (taken from Ericsson, Hendry and Prestwich 
    % Scandinavian Journal of Economics (1998, p.296)    
    %--------------------------------------------------------------------
    nsim = 50000;
    rm   = zeros( nsim,1 ) + 0.1;
    u    = zeros( nsim,1 );
    v    = 0.005*randn( nsim,1 );

    for i = 5:nsim

        rm(i) = rm(i-1) + 0.45*(rm(i-1) - rm(i-2)) - 0.10*(rm(i-2) ...
              - 2*rm(i-3) + rm(i-4))- 2.54*(rm(i-1) - 0.1)*rm(i-1)^2 + v(i);

    end


    drm = trimr( rm,1,0 ) - trimr( rm,0,1 );


    %**********************************************************************
    %***
    %***     Generate graph
    %***
    %**********************************************************************

    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(4);
    clf;

    hist( rm,51 );
    ylabel(' Frequency ');
    xlabel(' Midpoint ');
    %xlim([-5,5]);
    %set(gca,'XTick',-4:1:4);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    box off;
    axis tight;
    
    % Print the tex file to the relevant directory
    laprint(4,'multiple','options','factory');

    
end




