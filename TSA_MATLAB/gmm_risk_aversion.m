% ========================================================================
%
% GMM estimation of risk aversion using the Mehra-Prescott data
%    
% ========================================================================
function gmm_risk_aversion( )

    clear all;
    clc;
    
    % Load Mehra-Prescott data (annual for the period 1889 to 1978)
    load equity_mp.mat;
    dat = equity_mp;
                                          
    r  = dat(:,2);      %     Nominal risk free rate, yield percentage p.a.   
    p  = dat(:,3);      %     Price of consumption goods                      
    c  = dat(:,4);      %     Real per capita consumption                     
    d  = dat(:,5);      %     Real dividend                                   
    s  = dat(:,8);      %     Real stock price                                

    rs = (s(2:91) + d(1:90) - s(1:90))./s(1:90); 
    rb = (1 + r(1:90)./100).*((p(1:90))./p(2:91)) -1; 
    rc = (c(2:91) - c(1:90))./c(1:90);

    % Instruments
    const = ones(length(rs)-1,1);
    inst  = [ const   trimr(rc,0,1) ];
    %inst = [ const   trimr(rc,0,1)   trimr(rb,0,1) ];
    %inst = [ const   trimr(rc,0,1)   trimr(rb,0,1)   trimr(rs,0,1) ] ;
    
    % Adjust sample size
    rs = trimr(rs,1,0);
    rb = trimr(rb,1,0);
    rc = trimr(rc,1,0);
    t  = length(rs);

    % Identity weighting matrix
    flag        = 0;
    ops         = optimset('LargeScale','off','Display','off');
    theta0      = [ 0.9 ; 1.0 ];
    [theta,qmin] = fminunc(@(theta) q(theta,rc,rb,rs,inst,flag),theta0,ops);
    
    dof   = size(inst,2)-length(theta);
    Jstat = 2*t*qmin;

    disp('Results using identity weighting matrix ');
 	disp(['Risk aversion parameter estimate  = ', num2str(theta(2)) ]);
    disp(['J test statistic                  = ', num2str(Jstat) ]);
    
    if dof > 0.0; 
        disp(['p-value                       = ', num2str(1-chi2cdf(Jstat,dof)) ]);
    end

    
    % Optimal weighting matrix
    flag        = 1;
    [theta,qmin] = fminunc(@(theta) q(theta,rc,rb,rs,inst,flag),theta0,ops);

    dof   = size(inst,2)-length(theta);
    Jstat = 2*t*qmin;
    
    disp(' ');
    disp('Results using optimal weighting matrix ');
 	disp(['Risk aversion parameter estimate  = ', num2str(theta(2)) ]);
    disp(['J test statistic                  = ', num2str(Jstat) ]);
    
    if dof > 0.0; 
        disp(['p-value                       = ', num2str(1-chi2cdf(Jstat,dof)) ]);
    end

end
%
%------------------------- Functions -------------------------------------%
%
%-------------------------------------------------------------------------%
% GMM objective function   
%-------------------------------------------------------------------------%   
function ret = q(theta,rc,rb,rs,inst,flag)
        
    beta = theta(1);
    gam  = theta(2);
    m1   = bsxfun( @times,(beta*(1 + rc).^(-gam).*(1 + rb) - 1), inst );
    m2   = bsxfun( @times,(beta*(1 + rc).^(-gam).*(1 + rs) - 1), inst );
    m    = [ m1 m2 ];
    g    = mean(m)';
    
    if flag 
        w = m'*m/length(m);
    else
        w = eye(size(m,2));
    end
    ret  = 0.5*g'*inv(w)*g;
end