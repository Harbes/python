%==========================================================================
%
%   Simulation example to reproduce the asymptotic distribution of the 
%   MLE estimator for the regression model with autocorrelation.
%
%==========================================================================
function auto_distribution(  )

    clear all;
    clc;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234) )

    % Set parameter values
    beta0  = 1.0;                                           
    beta1  = 1.0;                                                     
    rho1   = 0.6;                                                       
    sig2   = 10;                 
    theta0 = [ beta0 ; beta1 ; rho1 ; sig2 ];  % Start values

    t       = 500;                                                                     
    ndraws  = 5000;
    options = optimset('LargeScale', 'off', 'Display', 'final');


    % Simulate a regression model with an AR(1) disturbance term 
    
    x = rand(t,1) - 0.5;  % x fixed in repeated samples
    
    %exact  = zeros(ndraws,1);
    cond   = zeros(ndraws,1);
    
    for k = 2:ndraws

        v = sqrt(sig2)*randn(t,1);   
        u = recserar( v , sqrt(1/(1-rho1^2))*v(1) , rho1 );             
        y = beta0 + beta1*x + u;
           
        % Exact MLE
        %[ tmp ]  =  fminunc(@(p) negloge(p,y,x),theta0,options);        
        %exact(k) = tmp(2);  
        
        % Conditional MLE
        
        [ tmp ] =  fminunc(@(p) neglogc(p,y,x),theta0,options);        
        cond(k) = tmp(2); 
     
    end
    
    mse       = mean( (cond - beta1).^2 );
    ssq       = sum( ( trimr(x,1,0) - rho1*trimr(x,0,1) ).^2 );                  
    i_beta    = ssq/sig2;             % Inforation matrix of beta hat    
    var_beta1 = inv(i_beta);          % Asymptotic variance              

    disp( [ 'Sample size                       = ' num2str(t) ] );
    disp( [ 'Sum of squares                    = ' num2str(ssq) ] );
    disp( [ 'Asymptotic variance (theoretical) = ' num2str(var_beta1) ] );
    disp( [ 'Asymptotic variance (simulated)   = ' num2str(mse) ] );
    
    %********************************************************************
    %***
    %***     Generate graph
    %***
    %********************************************************************

    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;
    
    z  = ( cond - beta1 )/sqrt( mse );
    tt = -5:0.1:5;
    
    [fcdf,xx] = ecdf(z); 
    [f,bins]  = ecdfhist(fcdf,xx,31); 
    bar(bins,f,'hist');
    
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k');
    
    box off
    hold on
    
    plot(tt,normpdf(tt),'-k','LineWidth',0.75)
    ylabel('$f(z)$');
    xlabel('$z$');
    set(gca,'YTick',[0.1 0.2 0.3 0.4 0.5])
    set(gca,'XTick',[-5 -4 -3 -2 -1 0 1 2 3 4 5])     
    
    hold off
    
    % Print the tex file to the relevant directory
    %laprint(1,'autodist','options','factory');

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
% Exact log-likelihood function
%-------------------------------------------------------------------------
function lf = negloge(b,y,x)

     beta0 = b(1);                                              
     beta1 = b(2);
     rho1  = tanh(b(3));    % Stay in the unit circle    
     sig2  = abs(b(4));     % Variance is positive
     fac   = 1 - rho1^2;

     u    = y - beta0 - beta1*x;
     v    = trimr(u,1,0) - rho1*trimr(u,0,1);
     tmp  = -0.5*log(2*pi)-0.5*log(sig2)+0.5*log(fac)-0.5*(u(1))^2/(sig2/fac); 
     tmp1 = - 0.5*log(2*pi) - 0.5*log(sig2) - 0.5*v.^2/sig2;
     
     lf = - sum( [ tmp; tmp1 ] );

end
%-------------------------------------------------------------------------
% Conditional log-likelihood function
%-------------------------------------------------------------------------
function lf = neglogc(b,y,x)

     beta0 = b(1);                                              
     beta1 = b(2);
     rho1  = tanh(b(3));    % Stay in the unit circle    
     sig2  = abs(b(4));     % Variance is positive
      
     u    = y - beta0 - beta1*x;
     v    = trimr(u,1,0) - rho1*trimr(u,0,1); 
     tmp  = - 0.5*log(2*pi) - 0.5*log(sig2) - 0.5*v.^2/sig2;
     
     lf = - sum( tmp );

end

