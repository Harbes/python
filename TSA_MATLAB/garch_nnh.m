%=========================================================================
%
%   Estimating ARCH-NNH Models by ML (Han and Park, J. Ects. (2008))
%
%=========================================================================
function garch_nnh( )

    clear all
    clc
    
    % Frequency
    freq = 'daily'; % 'daily','weekly','monthly'
    
    % Load data  
    load han_park;
    
    if strcmp(freq,'daily')
        spread = dailydata(:,3);
        rets   = dailydata(:,4);
    elseif strcmp(freq,'weekly')
        spread = weeklydata(:,3);
        rets   = weeklydata(:,4);
    elseif strcmp(freq,'monthly')
        spread = monthlydata(:,3);
        rets   = monthlydata(:,4);
    else
        disp('Wrong frequency....');
    end 
    
    % Define variables  
    y = trimr(rets,1,0);           
    w = trimr(spread,0,1);   % Lag spreads (not in percentage)  
    t = length(y);

    % Estimate the GARCH model
    ops           = optimset('LargeScale','off','Display','off');
    start         = [ 0.0006 0.001 0.04 0.9 ];
    [ theta,lf1 ] = fminsearch(@(b) neglog(b,y),start,ops);
    
   
    clear start
    lf1   = -lf1;
    theta = abs(theta);
    
    disp(' ');
    disp('GARCH(1,1) Results');
    disp(['mu                                      = ',num2str(theta(1)) ]);   
    disp(['alpha_0                                 = ',num2str(theta(2)) ]);
    disp(['alpha_1                                 = ',num2str(theta(3)) ]);
    disp(['beta_1                                  = ',num2str(theta(4)) ]);
    disp(['Log-likelihood function (unrestricted)  = ',num2str(lf1) ]);    

    % Estimate the GARCH-NNH model 
    flag = 'ARCH'; % 'ARCH' or 'GARCH'
    
    if strcmp(flag,'ARCH')
        start         = [ 0.05 0.2 0.1 1.0 ];
        [ theta1,lfn1 ] = fminunc(@(b) neglognnh(b,y,w,'ARCH',0),start,ops);
    
        clear start
        lfn1   = -lfn1;
        theta1 = abs(theta1); 
        disp(' ');
        disp('ARCH(1)-NNH Results');
        disp(['mu                                      = ',num2str(theta1(1)) ]);   
        disp(['alpha                                   = ',num2str(theta1(2)) ]);
        disp(['lambda                                  = ',num2str(theta1(3)) ]);
        disp(['phi                                     = ',num2str(theta1(4)) ]);
           
        start         = [ 0.05 0.2 0.1 ];
        [ ~,lfn0 ] = fminunc(@(b) neglognnh(b,y,w,'ARCH',1),start,ops);

        clear start
        lfn0 = -lfn0;
        
        % LR test
        disp(['Log-likelihood function (unrestricted)  = ',num2str(lfn1) ]); 
        disp(['Log-likelihood function (restricted)    = ',num2str(lfn0) ]); 

        lr = -2*t*(lfn0 - lfn1);
        disp(['LR test        = ',num2str(lr) ]);
        disp(['p-value        = ',num2str(1-chi2cdf(lr,1)) ]);

    elseif strcmp(flag,'GARCH')
        start         = [ 0.05 0.00 0.2 0.5 0.1 1.0 ];
        [ theta1,lfn1 ] = fminsearch(@(b) neglognnh(b,y,w,'GARCH',0),start,ops);
    
        clear start
        lfn1 = -lfn1;
        theta1 = abs(theta1);
        
        disp(' ');
        disp('GARCH(1)-NNH Results');
        disp(['mu                                      = ',num2str(theta1(1)) ]);   
        disp(['alpha_0                                 = ',num2str(theta1(2)) ]);
        disp(['alpha_1                                 = ',num2str(theta1(3)) ]);
        disp(['beta_1                                  = ',num2str(theta1(4)) ]);
        disp(['lambda                                  = ',num2str(theta1(5)) ]);
        disp(['phi                                     = ',num2str(theta1(6)) ]);
    
        start         = [ 0.05 0.00 0.2 0.5 0.1 ];
        [ ~,lfn0 ] = fminsearch(@(b) neglognnh(b,y,w,'GARCH',1),start,ops);

        
        disp(['Log-likelihood function (unrestricted)  = ',num2str(lfn1) ]); 
        disp(['Log-likelihood function (restricted)    = ',num2str(lfn0) ]); 

        lr = -2*t*(lfn0 - lfn1);
        disp(['LR test        = ',num2str(lr) ]);
        disp(['p-value        = ',num2str(1-chi2cdf(lr,1)) ]);


    end
        
    
end
%
%--------------------------- Functions  ----------------------------------
% 
%-------------------------------------------------------------------------
% Likelihood function for a GARCH(1,1) model
%-------------------------------------------------------------------------
function lf = neglog( b,y )
    
    b = abs(b);
    t = length( y );
    u = y-b(1); 
    h = std( y )^2*ones( t,1 );
    
    for i = 2:t
        
        h(i) = b(2) + b(3)*u(i-1)^2 + b(4)*h(i-1);

    end
    f  = - 0.5*log( 2*pi ) - 0.5*log( h ) - 0.5*(u./sqrt( h )).^2;
    lf = -mean( f );
    
end
%-------------------------------------------------------------------------
% Likelihood function for a GARCH-NNH(1,1) model
%-------------------------------------------------------------------------
function lf = neglognnh( b,y,w,flag,restrict )
    
    b = abs(b);  
    t = length( y );
    u = y-b(1); 
    h = std( y )^2*ones( t,1 );
    
    if strcmp(flag,'ARCH')      
        for i = 2:t      
            if restrict
                h(i) = b(2)*u(i-1)^2 + b(3)*abs(w(i));
            else
                h(i) = b(2)*u(i-1)^2 + b(3)*abs(w(i)).^b(4);
            end
        end
    elseif strcmp(flag,'GARCH')       
        for i = 2:t    
            if restrict
                h(i) = b(2)+b(3)*u(i-1)^2 + b(4)*h(i-1)+ b(5)*abs(w(i));
            else
                h(i) = b(2)+b(3)*u(i-1)^2 + b(4)*h(i-1)+ b(5)*abs(w(i)).^b(6);
            end
        end        
    else
        disp('Wrong flag in call to GARCH-NNH');
    end
    
    f  = - 0.5*log( 2*pi ) - 0.5*log( h ) - 0.5*(u./sqrt( h )).^2;
    lf = -mean( f );
    
end
