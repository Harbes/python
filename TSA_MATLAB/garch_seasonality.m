%=========================================================================
%
%   Estimate a GARCH(1,1) model of daily equity returns with seasonality 
%
%=========================================================================
function garch_seasonality(  )

    clear all;
    clc;
    
    % Load daily equity indices 
    % 5 January 1989 to 31 December 2007
    load equity

    equity  = ftse;   %  ftse,dow,nikkei        

    %     Compute the percentage return 
    y = 100*(trimr(log(equity),1,0) - trimr(log(equity),0,1));    
    y = y - mean(y);      
    t = length(y);

    % Trim dummy variables
    dtue = trimr(dum_tuesday,1,0);                         
    dwed = trimr(dum_wednesday,1,0); 
    dthu = trimr(dum_thursday,1,0); 
    dfri = trimr(dum_friday,1,0);                                      

    % Estimate the GARCH(1,1) model with seasonality     
    ops     = optimset( 'LargeScale','off','Display','off' );
    start   = [ 0.05 0.1 0.9 1 1 1 1 ];
    [~,lf1] = fminunc(@(b) neglog(b,y,dtue,dwed,dthu,dfri),start,ops);
    clear start
    
    lf1 = -lf1;

    disp(['Likelihood function (unrestricted)  = ',num2str(lf1) ]);

    % Estimate the GARCH(1,1) model without seasonality 
    start = [ 0.05 0.1 0.9 ];
    [~,lf0] = fminunc(@(b) neglog0(b,y),start,ops);
    
    lf0 = -lf0;

    disp(['Likelihood function (restricted)    = ',num2str(lf0) ]);

    % LR test
    lr  = -2*t*(lf0 - lf1);
    pv = 1-chi2cdf(lr,5);

	disp(' ');
    disp(['LR statistic           = ',num2str(lr) ]);
    disp(['p-value                = ',num2str(pv) ]);

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  GARCH(1,1) - unrestricted
%-------------------------------------------------------------------------
function lf = neglog(b,y,dtue,dwed,dthu,dfri)

    u = y;  
    d = b(4)^2*dtue + b(5)^2*dwed + b(6)^2*dthu + b(7)^2*dfri;                 
    h = recserar( b(1)^2 + b(2)^2*trimr([0.0;u.^2],0,1) + d,std(u)^2,b(3)^2);                     
    z = u./sqrt(h);                                                                       
    f =  - 0.5*log(2*pi) - 0.5*log(h) - 0.5*z.^2;
    
    lf = -mean( f );
end
%-------------------------------------------------------------------------
%  GARCH(1,1) - restricted
%-------------------------------------------------------------------------
function lf = neglog0(b,y)

    u = y;                                                                                         
    h = recserar(b(1)^2 + b(2)^2*trimr([0.0;u.^2],0,1),std(u)^2,b(3)^2);                    
    z = u./sqrt(h);                                                                          
    f =  - 0.5*log(2*pi) - 0.5*log(h) - 0.5*z.^2 ;
    
    lf = -mean( f );

end


