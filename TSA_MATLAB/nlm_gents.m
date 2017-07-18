%=========================================================================
%
%   Program to estimate the GENTS model using ftse equity returns
%
%=========================================================================
function nlm_gents( )

    clear all
    clc
    
    % Load data: daily returns on ftse 20 Nov. 1973 to 23 July 2001
    % Scale by 100 to convert to percentages   
    load ftse.mat              
    y = 100*data;
    
    %Treat outliers at observations 222, 233, 3527, 3528 and 3529
    dum1       = zeros(length(y),1);
    dum1(222)  = 1;
    dum2       = zeros(length(y),1);
    dum2(233)  = 1;
    dum3       = zeros(length(y),1);
    dum3(3527) = 1;
    dum4       = zeros(length(y),1);
    dum4(3528) = 1;
    dum5       = zeros(length(y),1);
    dum5(3529) = 1;

    % Rename y as the residuals from a dummy variable regression
    d = [ ones(length(y),1)   dum1   dum2   dum3   dum4   dum5 ];
    y = y - d*(d\y);                                         

    % Current and lagged returns
    x = trimr(y,0,1);                                                 
    y = trimr(y,1,0);                              
    t = length(y);

    % Estimate model
    ops     = optimset('LargeScale','off','Display','iter');
    %theta_0 = [10 ; 1 ; 1 ];
    
    theta_0 = [      0.928834745978081
                    -0.700513899637738
                     0.695374988151895 ];

    %theta_0 = [14 -1 0.1];             
    
    [theta,~,~,~,~,hess] = fminunc(@(b) neglog(b,y,x),theta_0,ops);
    % [theta,~ ] = fminsearch(@(b) neglog(b,y,x),theta_0,ops);
   
    gam = theta(1);
    vc  = (1/t)*inv(hess);

    % Perform Wald test of skewness   
    r = [  0   1   0 ;
           0   0   1  ]; 
    q = [ 0 ; 0 ];
    w = (r*theta - q)'*inv(r*vc*r')*(r*theta - q);

    disp(' ' );
    disp(['Wald test of skewness = ',num2str(w) ]);
    disp(['p-value               = ',num2str(1-chi2cdf(w,2))]);

    % Compute conditional moments    
    y_mean = zeros(t,1);
    y_var  = zeros(t,1);
    y_skew = zeros(t,1);
    y_kurt = zeros(t,1);
    
    gam = theta(1);
    t1  = theta(2)*x;                      
    t2  = -0.5*(1+gam^2)*ones(t,1);
    t3  = theta(3)*x;
    t4  = -0.5*ones(t,1);
    
    [~,eta] = neglog(theta,y,x);


    for i = 1:t
        y_mean(i) = quadgk(@(s) gst1(s,gam,t1(i),t2(i),t3(i),t4(i),eta(i)),-15,15);
        y_var(i)  = quadgk(@(s) gst2(s,gam,y_mean(i),t1(i),t2(i),t3(i),t4(i),eta(i)),-15,15);
%       y_skew(i) = quadgk(@(s) gst3(s,gam,y_mean(i),t1(i),t2(i),t3(i),t4(i),eta(i)),-15,15);
%       y_kurt(i) = quadgk(@(s) gst4(s,gam,y_mean(i),t1(i),t2(i),t3(i),t4(i),eta(i)),-15,15);
    end

    figure(1)
    plot(1:t,y_mean); 
    figure(2)
    plot(1:t,y_var);
%     figure(3)
%     plot(1:t,y_skew./y_var.^1.5);
%     figure(4)
%     plot(1:t,y_kurt./y_var.^2);
end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Log-likelihood function
%-------------------------------------------------------------------------
function [lf,eta] = neglog(b,y,x)

    t   = length(y);
    gam = b(1);                          
    t1  = b(2)*x;                      
    t2  = -0.5*(1+gam^2)*ones(t,1);
    t3  = b(3)*x;
    t4  = -0.5*ones(t,1);
    
    eta = zeros(t,1);
    lnl = zeros(t,1);
    % Range of integration includes the max and min values of yt  
    for i = 1:t      
        eta(i)  = log( quadgk(@(s) gst(s,gam,t1(i),t2(i),t3(i),t4(i)),-15,15,'RelTol',1e-8,'AbsTol',1e-12) );  
        lnl(i) = t1(i)*atan(y(i)/gam) + t2(i)*log(gam^2+y(i)^2) + t3(i)*y(i) + t4(i)*y(i)^2 - eta(i);

    end
    lf = -mean( lnl );

end
%-------------------------------------------------------------------------
%  Compute eta for the generalised Student t distribution
%-------------------------------------------------------------------------
function f = gst(s,gam,t1,t2,t3,t4)

     f = exp( t1*atan(s/gam) + t2*log(gam^2 + s.^2) + t3*s + t4*s.^2 );
    
end
%-------------------------------------------------------------------------
%  Compute the conditional mean using the generalised Student t distribution
%-------------------------------------------------------------------------
function f = gst1(s,gam,t1,t2,t3,t4,eta)

     f = s.*exp( t1*atan(s/gam) + t2*log(gam^2 + s.^2) + t3*s + t4*s.^2  - eta);
    
end
%-------------------------------------------------------------------------
%  Compute the conditional variance using the generalised Student t distribution
%-------------------------------------------------------------------------
function f = gst2(s,gam,mm,t1,t2,t3,t4,eta)

    tmp = (s - mm).^2;
    f   = tmp.*exp( t1.*atan(s/gam) + t2.*log(gam^2 + s.^2) + t3*s + t4*s.^2  - eta);
    
end
%-------------------------------------------------------------------------
%  Compute the conditional skewness using the generalised Student t distribution
%-------------------------------------------------------------------------
function f = gst3(s,gam,mm,t1,t2,t3,t4,eta)

    tmp = (s - mm).^3;
    f   = tmp.*exp( t1.*atan(s/gam) + t2*log(gam^2 + s.^2) + t3*s + t4*s.^2 - eta);
    
end
%-------------------------------------------------------------------------
%  Compute the conditional kurtosis using the generalised Student t distribution
%-------------------------------------------------------------------------
function f = gst4(s,gam,mm,t1,t2,t3,t4,eta)

    tmp = (s - mm).^4;
    f   = tmp.*exp( t1*atan(s/gam) + t2*log(gam^2 + s.^2) + t3*s + t4*s.^2 - eta);
    
end
