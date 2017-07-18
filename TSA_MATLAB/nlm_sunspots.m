%=========================================================================
%
%   Sunspot data
%
%=========================================================================
function nlm_sunspots( )

    clear all
    clc
    format short
    
    load sunspot.mat
    
    y = sunspot(:,1);
    t = length(y);

    figure(1)
    clf;
    plot(1:length(y),y,'-k');
    ylabel('Average Monthly Sunspots');
    axis tight
    box off

    % Perform LST linearity test   
    [lm,dof] = tar_test(y,1);          
    disp(['LM statistic (test 1) = ', num2str(lm) ]);
    disp(['Degrees of freedom    = ', num2str(dof) ]);
    disp(['p-value               = ', num2str(1-chi2cdf(lm,dof))]);
    disp(' ');

    [lm,dof] = tar_test(y,2);          
    disp(['LM statistic (test 2) = ', num2str(lm) ]);
    disp(['Degrees of freedom    = ', num2str(dof) ]);
    disp(['p-value               = ', num2str(1-chi2cdf(lm,dof))]);
    disp(' ');

    [lm,dof] = tar_test(y,3);          
    disp(['LM statistic (test 3) = ', num2str(lm) ]);
    disp(['Degrees of freedom    = ', num2str(dof) ]);
    disp(['p-value               = ', num2str(1-chi2cdf(lm,dof))]);
    disp(' ');

    [lm,dof] = tar_test(y,4);          
    disp(['LM statistic (test 4) = ', num2str(lm) ]);
    disp(['Degrees of freedom    = ', num2str(dof) ]);
    disp(['p-value               = ', num2str(1-chi2cdf(lm,dof))]);
    disp(' ');


    % Parameters
    d   = 2;        % Delay parameter     
    gam = 1;        % Adjustment parameter   1, 5, 10, 50, 100       

    % Estimate the model
    ops     = optimset('LargeScale','off','Display','off' );   

    theta_0 = [ ones(11,1)*0.1 ; 10 ];

    [theta,lnl,~,~,~,hess] = fminunc(@(b) neglog(b,y,d,gam),theta_0,ops);

    vc = (1/t)*inv(hess);
    
    disp(' ')
    disp(['Log-likelihood    = ',num2str(-lnl);])
    disp(' Parameters    Std. Errors')
    disp( [theta sqrt(diag(vc))] )

end

%-------------------------------------------------------------------------
%  Compute the LM statistic to test a threshold autoregressive model 
%  assuming one lag in the auxiliary model
%-------------------------------------------------------------------------
function [ lm,dof ] = tar_test(yvar,p)

    % First stage regression      
    y = trimr(yvar,1,0);
    x = [ones(size(y,1),1) , trimr(yvar,0,1)];                     
    k = size(x,2);
    u = y - x*(x\y);                                        

    % Second stage regression      
    if p == 1
        x = [x  x(:,2).*(x(:,2))]; 
    elseif p == 2
        x = [x  x(:,2).*(x(:,2)).^1  x(:,2).*(x(:,2).^2)]; 
    elseif p == 3
        x = [x , x(:,2).*(x(:,2).^1) , x(:,2).*x(:,2).^2 , x(:,2).*x(:,2).^3]; 
    elseif p == 4
        x = [x , x(:,2).*(x(:,2).^1) , x(:,2).*(x(:,2).^3)]; 
    end     
    e = u - x*(x\u);                                        

    % Compute LM statistic        
    r2  = 1 - sum(e.^2)/sum( (u-mean(u)).^2 );
    lm  = length(y)*r2;
    dof = size(x,2) - k;  
    
end
%-------------------------------------------------------------------------
%  Negative log-likelihood function
%-------------------------------------------------------------------------
function lf = neglog(b,y,d,gam)

     c  = b(12);
     w  =  1./( 1 + exp(-gam*(trimr(y,6-d,d) - c)));
     
     % List of variables in regime 1       
     x1 = [ ones(length(y)-6,1)   trimr(y,5,1)   trimr(y,4,2)   trimr(y,3,3)   trimr(y,2,4)   trimr(y,1,5)   trimr(y,0,6) ];    

     % List of variables in regime 2      
     x2 = [ ones(length(y)-6,1)   trimr(y,5,1)   trimr(y,4,2)   trimr(y,3,3) ];                                             

     u  = trimr(y,6,0) - x1*b(1:7) - x2*b(8:11).*w;

     sig2 = u'*u/length(u);             
     lnl  = - 0.5*log(2*pi) - 0.5*log(sig2) - 0.5*u.^2/sig2;                                                 
     lf   = -mean( lnl );

end
