%=========================================================================
%
%   Program to estimate a bivariate threshold model
%   Anderson, Anthanasopoulos Vahid (2007, 22, 63–87).
%
%   Uses quarterly data on percentage growth rate in real gdp and 
%   interest rate spread (10yr - 3mth). Source: JAE Archive.
%
%=========================================================================
function nlm_g7(  )

    % Set country
    % Canada, France, Germany, Italy, Japan, UK, US
    country = 'US';    
  
     % Get the data 
     y    = getData( country );
     y1   = y(:,1);                 % Growth rate of GDP
     y2   = y(:,2);                 % Spread
     t    = length( y1 );

    % Perform bivariate LM nonlinearity tests on GDP 
    dep          = 1;
    trans        = 2;
    [ lm1,dof1 ] = btar_test( y,dep,trans,1 );
    [ lm2,dof2 ] = btar_test( y,dep,trans,2 );
    [ lm3,dof3 ] = btar_test( y,dep,trans,3 );
    [ lm4,dof4 ] = btar_test( y,dep,trans,4 );

    disp('LM tests of nonlinearity in GDP');
    disp('-------------------------------');
    disp(['LM statistic (Model 1) =  ', num2str(lm1) ]);
    disp(['p-value                =  ', num2str(1 - chi2cdf(lm1,dof1)) ]);
    disp(' ') 
    disp(['LM statistic (Model 2) =  ', num2str(lm2) ]);
    disp(['p-value                =  ', num2str(1 - chi2cdf(lm2,dof2)) ]);
    disp(' ') 
    disp(['LM statistic (Model 3) =  ', num2str(lm3) ]);
    disp(['p-value                =  ', num2str(1 - chi2cdf(lm3,dof3)) ]);
    disp(' ') 
    disp(['LM statistic (Model 4) =  ', num2str(lm4) ]);
    disp(['p-value                =  ', num2str(1 - chi2cdf(lm4,dof4)) ]);
    disp(' ') 


    % Perform bivariate LM nonlinearity tests on the Spread
    dep          = 2;
    trans        = 2;
    [ lm1,dof1 ] = btar_test( y,dep,trans,1 );
    [ lm2,dof2 ] = btar_test( y,dep,trans,2 );
    [ lm3,dof3 ] = btar_test( y,dep,trans,3 );
    [ lm4,dof4 ] = btar_test( y,dep,trans,4 );

    disp('LM tests of nonlinearity in the Spread');
    disp('--------------------------------------');
    disp(['LM statistic (Model 1) =  ', num2str(lm1) ]);
    disp(['p-value                =  ', num2str(1 - chi2cdf(lm1,dof1)) ]);
    disp(' ') 
    disp(['LM statistic (Model 2) =  ', num2str(lm2) ]);
    disp(['p-value                =  ', num2str(1 - chi2cdf(lm2,dof2)) ]);
    disp(' ') 
    disp(['LM statistic (Model 3) =  ', num2str(lm3) ]);
    disp(['p-value                =  ', num2str(1 - chi2cdf(lm3,dof3)) ]);
    disp(' ') 
    disp(['LM statistic (Model 4) =  ', num2str(lm4) ]);
    disp(['p-value                =  ', num2str(1 - chi2cdf(lm4,dof4)) ]);
    disp(' ') 


    % Estimate the bivariate model using BFGS  
    final  = zeros(22,7);
    maxlag = 2; % Maximum number of lags in bivariate model 
    start  = 0.1*ones(22,1);   %     Starting parameter estimates        
    ops    = optimset('LargeScale','off','Display','iter', ...
                       'MaxIter',200,'MaxFunEvals',5000');

    [ bhat,lf,~,~,~,hess ] = fminunc( @(b) neglog(b,t,y1,y2,maxlag),start,ops);
 
    lf    = -lf;
    omega = (1/t)*inv(hess);

    disp(' ') 
    disp( country )
    disp([ 'Log-likelihood function= ', num2str(lf) ]);
    disp('Covariance matrix')
    disp(omega);
    
    % For plotting and simulation use GAUSS estimates
    final = xlsread('FinalEstimates.xlsx','A1:G22');

   % Simulate all models 
    ystart = [1 ; 1];   % Start values      
 
    % Canada
    tmp = btarsim(final(:,1),ystart,10000,maxlag);
    y1s_canada = tmp(:,1);
    y2s_canada = tmp(:,2);
    
    % France
    tmp = btarsim(final(:,2),ystart,10000,maxlag);
    y1s_france = tmp(:,1);
    y2s_france = tmp(:,2);
   
    % Germany
    tmp = btarsim(final(:,3),ystart,10000,maxlag);
    y1s_germany = tmp(:,1);
    y2s_germany = tmp(:,2);

    % Italy
    tmp = btarsim(final(:,4),ystart,10000,maxlag);
    y1s_italy = tmp(:,1);
    y2s_italy = tmp(:,2);
    
    % Japan
    tmp = btarsim(final(:,5),ystart,10000,maxlag);
    y1s_japan = tmp(:,1);
    y2s_japan = tmp(:,2);

    % UK
    tmp = btarsim(final(:,6),ystart,10000,maxlag);
    y1s_uk = tmp(:,1);
    y2s_uk = tmp(:,2);
    
    % US
    tmp = btarsim(final(:,7),ystart,10000,maxlag);
    y1s_us = tmp(:,1);
    y2s_us = tmp(:,2);
    
%**********************************************************************
%***
%***     Generate graph
%***
%**********************************************************************
subplot(7,3,1)
    xlabel( 't ');
    ylabel('Canada-output ');
    plot(y1s_canada(1:100));
subplot(7,3,2)   
    xlabel( 't ');
    ylabel('Canada-spread ');
    plot(y2s_canada(1:100));
subplot(7,3,3)  
    xlabel('Canada-output');
    ylabel('Canada-spread');
    plot(y1s_canada,y2s_canada);
subplot(7,3,4) 
    xlabel('t ');
    ylabel('France-output ');
    plot(y1s_france(1:100));
subplot(7,3,5)
    xlabel('t ');
    ylabel('France-spread ');
    plot(y2s_france(1:100));
subplot(7,3,6)
    xlabel('France-output ');
    ylabel('France-spread ');
    plot(y1s_france,y2s_france);
subplot(7,3,7)
    xlabel(' t ');
    ylabel('Germany-output ');
    plot(y1s_germany(1:100));
subplot(7,3,8)
    xlabel(' t ');
    ylabel('Germany-spread ');
    plot(y2s_germany(1:100));
subplot(7,3,9)
    xlabel('Germany-output ');
    ylabel('Germany-spread ');
    plot(y1s_germany,y2s_germany);
subplot(7,3,10)
    xlabel(' t ');
    ylabel('Italy-output ');
    plot(y1s_italy(1:100));
subplot(7,3,11)    
    xlabel(' t ');
    ylabel('Italy-spread ');
    plot(y2s_italy(1:100));
subplot(7,3,12)  
    xlabel('Italy-output ');
    ylabel('Italy-spread ');
    plot(y1s_italy,y2s_italy);
subplot(7,3,13)   
    xlabel(' t ');
    ylabel('Japan-output ');
    plot(y1s_japan(1:100));
subplot(7,3,14)   
    xlabel(' t ');
    ylabel('Japan-spread ');
    plot(y2s_japan(1:100));
subplot(7,3,15)
    xlabel('Japan-output ');
    ylabel('Japan-spread ');
    plot(y1s_japan,y2s_japan);
subplot(7,3,16)
    xlabel(' t ');
    ylabel('UK-output ');
    plot(y1s_uk(1:100));
subplot(7,3,17)
    xlabel(' t ');
    ylabel('UK-apread ');
    plot(y2s_uk(1:100));
subplot(7,3,18)
    xlabel('UK-output ');
    ylabel('UK-spread ');
    plot(y1s_uk,y2s_uk);
subplot(7,3,19)
    xlabel('t ');
    ylabel('US-output ');
    plot(y1s_us(1:100));
subplot(7,3,20)
    xlabel(' t ');
    ylabel('US-spread ');
    plot(y2s_us(1:100));
subplot(7,3,21)
    xlabel('US-output ');
    ylabel('US-spread ');
    plot(y1s_us,y2s_us);
    
end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
% Unrestricted likelihood function 
%-------------------------------------------------------------------------
function lf = neglog(b,nobs,y1,y2,maxlag)

     v1 = zeros(nobs,1);     
     v2 = zeros(nobs,1);
     lt = zeros(nobs,1);

     for t = maxlag+1:nobs

        w1 = 1/( 1 + exp(-100*(y1(t-2) - b(11))  ) );        
        v1(t) = y1(t) - ( b(1)  + b(2)*y1(t-1)  + b(3)*y1(t-2)  + b(4)*y2(t-1)  + b(5)*y2(t-2) ) ...
                      - ( b(6) + b(7)*y1(t-1) + b(8)*y1(t-2) + b(9)*y2(t-1) + b(10)*y2(t-2) )*w1;

        w2 = 1/( 1 + exp(-100*(y1(t-1) - b(22)) ) );         
        v2(t) = y2(t) - ( b(12) + b(13)*y1(t-1) + b(14)*y1(t-2) + b(15)*y2(t-1) + b(16)*y2(t-2) ) ...
                      - ( b(17) + b(18)*y1(t-1) + b(19)*y1(t-2) + b(20)*y2(t-1) + b(21)*y2(t-2) )*w2;
         
     end
     
     %     Exclude the first maxlag observations  
     v     = trimr([v1,v2],maxlag,0);        
     [n,k] = size(v);
     
     vc = v'*v/n;

    for t = 1:nobs-maxlag

      lt(t) = - k*0.5*log(2*pi) - 0.5*log(det(vc)) - 0.5*v(t,:)*inv(vc)*v(t,:)';

    end

    %  Log-likelihood is not defined for the last maxlag observations  
    lf = - mean( trimr(lt,0,maxlag) );      
  
end
%------------------------------------------------------------------------- 
%   Gets the data for a specfic country
%------------------------------------------------------------------------- 
function y = getData( flag )

    if strcmp(flag,'Canada')       % Canada
    
        y = xlsread('G7Data.xlsx','B8:C162');
    
    elseif strcmp(flag,'France')   % France    
    
        y = xlsread('G7Data.xlsx','D43:E157');
    
    elseif strcmp(flag,'Germany')    % Germany
    
        y = xlsread('G7Data.xlsx','F4:G162');
    
    elseif strcmp(flag,'Italy')    % Italy
    
        y = xlsread('G7Data.xlsx','H48:I162');
    
    elseif strcmp(flag,'Japan')    % Japan
    
        y = xlsread('G7Data.xlsx','J43:K160');
    
    elseif strcmp(flag,'UK')    % UK
    
        y = xlsread('G7Data.xlsx','L4:M162');
    
    else
    
        y = xlsread('G7Data.xlsx','N4:O162');
    
    end
    
end
%------------------------------------------------------------------------- 
%   Compute the LM statistic to test a bivariate TAR
%   model assuming one lag in the auxiliary model
%------------------------------------------------------------------------- 
function [ lm,dof ] = btar_test(yvar,dep,trans,itype)

    % itype = 1: regress yt on [constant ylag, ylag*ylag]       
    % itype = 2: regress yt on [constant ylag, ylag*ylag, ylag*ylag^2]   
    % itype = 3: regress yt on [constant ylag, ylag*ylag, ylag*ylag^2, ylag*ylag^3]  
    % itype = 4: regress yt on [constant ylag, ylag*ylag, ylag*ylag^3]  

    % First stage regression     
    y = trimr( yvar(:,dep),1,0);           % Choose the dependent variable 
    x = [ ones( length(y),1 ) trimr( yvar,0,1 ) ];                     
    k = size( x,2 );
    u = y - x*(x\y);                                       

    % Second stage regression    
    z = trimr(yvar(:,trans),0,1);           %     Choose the transition variable (lagged one period)

    if itype == 1;      
        
        tp = bsxfun(@times,x(:,2:k),z.^1);
        x   = [x , tp]; 
    
    elseif itype == 2; 
        
        tp = bsxfun(@times,x(:,2:k),z.^1);
        tt = bsxfun(@times,x(:,2:k),z.^2);
        x   = [x tp tt]; 
    
    elseif itype == 3; 
        
        tp = bsxfun(@times,x(:,2:k),z.^1);
        tt = bsxfun(@times,x(:,2:k),z.^2);
        tr = bsxfun(@times,x(:,2:k),z.^3);
        x   = [x tp tt tr]; 
    
    elseif itype == 4; 
        
        tp = bsxfun(@times,x(:,2:k),z.^1);
        tt = bsxfun(@times,x(:,2:k),z.^3);
        x   = [x tp tt]; 
    
    end
    e = u - x*(x\u);                                        

    % Compute LM statistic        
    t  = size(y,1);
    r2 = 1 - sum(e.^2)/sum( (u-mean(u)).^2 );
    lm = t*r2;

    dof = size(x,2) - k;  
  
end
%------------------------------------------------------------------------- 
%   Simulate the model 
%------------------------------------------------------------------------- 
function y = btarsim(b,ystart,nobs,maxlag)

     y1 = zeros(nobs,1)+ystart(1);     
     y2 = zeros(nobs,1)+ystart(2);

     for t = maxlag+1:nobs

        w1 = 1/( 1 + exp(-100*(y1(t-2) - b(11))  ) );        
        y1(t) =  ( b(1)  + b(2)*y1(t-1)  + b(3)*y1(t-2)  + b(4)*y2(t-1)  + b(5)*y2(t-2) ) ...
                      - ( b(6) + b(7)*y1(t-1) + b(8)*y1(t-2) + b(9)*y2(t-1) + b(10)*y2(t-2) )*w1;

        w2 = 1/( 1 + exp(-100*(y1(t-1) - b(22)) ) );         
        y2(t) =  ( b(12) + b(13)*y1(t-1) + b(14)*y1(t-2) + b(15)*y2(t-1) + b(16)*y2(t-2) ) ...
                      - ( b(17) + b(18)*y1(t-1) + b(19)*y1(t-2) + b(20)*y2(t-1) + b(21)*y2(t-2) )*w2;
         
     end
     
     %     Exclude the first maxlag observations  
     y     = trimr([y1 y2],maxlag,0);        

end
