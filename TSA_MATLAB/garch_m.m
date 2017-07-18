%=========================================================================
%
%   Estimating GARCH-M Models
%
%=========================================================================
function garch_m( )

    clear all
    clc
    
    % Load data from file
    load yields_us
    
    data = rdata;
    r3m = data(:,1);         % Independent variable is the 3 month yield
    r10 = data(:,6);         % Dependent variable is the 10 year yield 
    t   = length(r10);
    
    % Estimating the GARCH-M(1,1) model
    ops = optimset( 'LargeScale','off','Display','iter');  

    % Starting values are close to the maximum likelihood estimates         
    start = [   0.0489928275220468 
                0.9617406702015530 
                0.2753291619382821 
                2.2590829301737938 
                0.7763991081157382 
                0.1167781704457386 
                1.0354494753290766];
    
    [theta,lf1,~,~,~,hess] = fminunc( @(b) neglog(b,r3m,r10),start,ops);    

    clear start
    lf1 = -lf1;
    vc = (1/t)*inv(hess);
    
    disp(' ');
    disp('GARCH-M(1,1) Results');
    disp('Mean parameters')
    disp(['gamma0                                  = ',num2str(theta(4)) ]); 
    disp(['gamma1                                  = ',num2str(theta(5)) ]);
    disp(['varphi                                  = ',num2str(theta(6)) ]);
    disp(['rho                                     = ',num2str(theta(7)) ]);
    disp(['alpha_0                                 = ',num2str(theta(1)) ]);
    disp(['alpha_1                                 = ',num2str(theta(2)) ]);
    disp(['beta_1                                  = ',num2str(theta(3)) ]);
    disp(['Log-likelihood function (unrestricted)  = ',num2str(lf1) ]);    
    
    % Wald test for risk neutrality    
    w = (theta(6)-0.0)^2/vc(6,6);
    
    disp(['Wald test of varphi = 0    = ',num2str(w) ]); 
    disp(['p-value                    = ',num2str(1-chi2cdf(w,1)) ]);    
    
    % Wald test for risk-return relationship   
    w = (theta(7)-1.0)^2/vc(7,7);
    disp(['Wald test of rho = 1       = ',num2str(w) ]); 
    disp(['p-value                    = ',num2str(1-chi2cdf(w,1)) ]);    
  
end

%--------------------------- Subroutines ----------------------------------
% 
%-------------------------------------------------------------------------
% Log-likelihood function 
%-------------------------------------------------------------------------

function lf = neglog(b,r3m,r10)
   	 
    b = abs(b); 
    t = length( r3m );
    m  = zeros(t,1);
    u = zeros(t,1); 
    h = std( r10 - r3m )^2*ones( t,1 );
    
    for i = 2:t       
        h(i) = b(1) + b(2)*u(i-1)^2 + b(3)*h(i-1);
        m(i) = b(4) + b(5)*r3m(i) + b(6)*(h(i)^(0.5))^b(7);
        u(i) = r10(i) - m(i);
    end   
    z  = u./sqrt(h);     
	f  = -0.5*log(2*pi) - 0.5*log(h) - 0.5*z.^2;
    lf = -mean( f );
end
