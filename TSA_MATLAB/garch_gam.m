%=========================================================================
%
%  Efficiency of QMLE for a garch model with conditional gamma
%
%=========================================================================
function garch_gam( )

    clear all
    clc
    
    % Set random number generator
    RandStream.setGlobalStream( RandStream('mt19937ar','seed',12345) );

    t = 1000000;
    z    = randn(t,1);                       
    v    = z;                                                
    y    = z;
    s2   = ones(t,1);

    % Parameters
    a1 = 0.1;
    b1 = 0.8;
    c  = 50;

    %  Simulate the dgp of a garch model with conditional gamm disturbance     
    for i = 2:t

        z(i) = (gamrnd(c,1) - c)/sqrt(c);                             
        s2(i) = 1-a1-b1 + a1*v(i-1)^2 + b1*s2(i-1);         
        v(i) = z(i)*sqrt(s2(i));                            
        y(i) = v(i);                                        
    end
    
    % Analaytical expressions (based on Engle and Gonzalez-Rivera (JBES, 1991)     
    ds2a = recserar( [0 ; trimr(y.^2,0,1) - 1], 0.0 , b1 );
    ds2b = recserar( [0 ; trimr(s2,0,1)   - 1], 0.0 , b1 );
    ds2  = [ ds2a   ds2b ];

    tmp = bsxfun(@rdivide,ds2,s2);
    cov_h    = inv( 0.5*(tmp)'*(tmp)/t );
    beta2    = 3 + 6/c;
    cov_j    = inv( (beta2-1)*0.25*(tmp)'*(tmp)/t );
    cov_qmle = cov_h*inv(cov_j)*cov_h;  
    
    clear tmp

    tmp = bsxfun(@times,ds2,(c - c*z.^2)./( c + sqrt(c)*z ));
    dl  = bsxfun(@rdivide,-0.5*tmp,s2);
    j0    = dl'*dl/t;
    cov_0 = inv(j0);                                        

    % Compute relative efficiency of qmle (based on analytical derivatives)     

    disp('Relative efficiency of qmle using analytical derivatives')
    disp('    alpha1    beta1')
    disp( (diag(cov_0)./diag(cov_qmle))' ) 
    disp([' alpha1 = ',num2str(a1) ]);
    disp([' beta1  = ',num2str(b1) ]);
    disp([' shape  = ',num2str(c) ]);
    disp( ' ' )

    % GARCH(1,1) gamma distribution:  
    ops   = optimset( 'LargeScale','off','Display','off');  
    start = [a1 b1 c];
    [theta,~,~,~,~,hess] = fminunc( @(b) lnl(b,y),start,ops);
     
    g   = numgrad( @lnlt,theta',y );
    %vcg = g'*g/t;
    vch = inv(hess);
    tt  = diag(vch);

    % GARCH(1,1) normal distribution: qmle covariance
    start  = [a1 b1];
    [theta0,~,~,~,~,h0] = fminunc( @(b) lnl1(b,y),start,ops);
        
    ih   = inv(h0);
    g0   = numgrad( @lnlt1,theta0',y );
    j0   = g0'*g0/t;
    vc0  = ih*j0*ih;
    tt0  = diag(vc0);
        
    % Efficiency ratio based on numerical derivatives
    disp('Numerical results');
    disp( (tt(1:2)./tt0)' );
    
end
%
%--------------------------- Functions ----------------------------------
% 
%-------------------------------------------------------------------------
% Likelihood wrapper function
% This function calls the lnlt function and returns the average log-likelihood.
%-------------------------------------------------------------------------

function logl = lnl( b,y )

    logl = -mean( lnlt( b,y ) ); 
    
end

%-------------------------------------------------------------------------
% Likelihood function for a GARCH-gam(1,1) model
%-------------------------------------------------------------------------

function loglt = lnlt( b,y )
    
    v  = y;                                                                     
    s2 = recserar(1-b(1)-b(2) + b(1)*trimr([0.0; v.^2],0,1),std(y)^2,b(2));     
    c  = abs(b(3));                                             
    z  = v./sqrt(s2);                                                            
             
    loglt = 0.5*log(c) - log(gamma(c)) + (c-1)*log(sqrt(c)*z + c) ...
           - 0.5*log(s2) - sqrt(c)*z - c;
end

%-------------------------------------------------------------------------
% Likelihood wrapper function
% This function calls the lnlt function and returns the average log-likelihood.
%-------------------------------------------------------------------------

function logl = lnl1( b,y )

    logl = -mean( lnlt1( b,y ) ); 
    
end

%-------------------------------------------------------------------------
% Likelihood function for a GARCH-N(1,1) model
%-------------------------------------------------------------------------

function loglt1 = lnlt1( b,y )
    
    b = abs(b);
    t = length( y );
    u = y;  
    h = std( y )^2*ones( t,1 );
    
    for i = 2:t
        
        h(i) = 1-b(1)-b(2) + b(1)*u(i-1)^2 + b(2)*h(i-1);

    end
    loglt1 = - 0.5*log( 2*pi ) - 0.5*log( h ) - 0.5*(u./sqrt( h )).^2;
    
end

% 
% /**     True model: Log of the likelihood for a GARCH(1,1) model with conditional Student t disturbance     **/
% 
%     proc lnlt(b,y);
% 
%         local v,s2,c,z,const,lnl;
% 
% 
%     endp;
% 
% 
% 
% /**     Misspecified model: Log of the likelihood for a GARCH(1,1) model with conditional normal disturbance     **/
% 
%     proc lnlt_normal(b,y);
% 
%         local v,s2,z,lnl;
% 
%         v  = y;                                                                     /**     Disturbance                                     **/
% 
%         s2 = recserar(1-b[1]-b[2] + b[1]*trimr(0.0|v.^2,0,1),stdc(y)^2,b[2]);       /**     Variance                                        **/
% 
%         z = v./sqrt(s2);                                                            /**     Standardised disturbance                        **/
% 
%         lnl = -0.5*ln(2*pi) - 0.5*ln(s2) - 0.5*z.^2;                                /**     Log likelihood                                  **/
% 
%         retp( lnl );
% 
%     endp;
% 
% 
% 
% output off;