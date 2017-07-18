%=========================================================================
%
%   Simulating GARCH - Student t Model
%  
%=========================================================================
function garch_studt( )

    clear all
    clc
    
    % Set random number generator
    RandStream.setGlobalStream( RandStream('mt19937ar','seed',1234) )
    
    t  = 100000;
    
    z  = randn( t,1 );        
    u  = z;                         
    y  = z;
    h  = ones( t,1 );

    a1 = 0.1;
    b1 = 0.8;
    nu = [5,8,12]; 
    
    % Make only one call to tinv and rand
    w  = tinv(rand(t,3),ones(t,3)*diag(nu));
    
    % Parameter values first model
    for k = 1:length(nu)
    
        disp(' ');
        disp(['Degrees of freedom: nu  = ',num2str(nu(k)) ]);
        % Generate data 
        for i = 2:t
        
            z(i) = sqrt((nu(k)-2)/nu(k))*w(i,k);                % Student t disturbance
            h(i) = 1-a1-b1 + a1*u(i-1)^2 + b1*h(i-1);         % Conditional variance  
            u(i) = z(i)*sqrt( h(i) );                         % Disturbance term     
            y(i) = u(i);                                      % Returns
        end
  
        % Efficiency based on analytical derivatives Engle and Gonzalez-Rivera (JBES, 1991)
        ds2a = recserar( [0 ; trimr(y.^2,0,1)] - 1, 0.0 , b1 );
        ds2b = recserar( [ 0 ; trimr(h,0,1)]    - 1, 0.0 , b1 );
        ds2  = [ ds2a   ds2b ];
    
        tmp   = bsxfun(@rdivide,ds2,h); 
        tmp1  = tmp'*tmp/t;
        cov_h = inv( 0.5*tmp1 );
        beta2 = 3 + 6/(nu(k)-4);
        cov_j = inv( (beta2-1)*0.25*tmp1 );
        
        % Analytical qmle covariance matrix of the misspecified model 
        cov_qmle = cov_h*inv(cov_j)*cov_h;                      

        % Analytical covariance matrix of the true model
        tmp2  = bsxfun(@times,ds2,(1 - ((nu(k)+1)*y.^2)./( y.^2 + h*(nu(k)-2) ) ));
        dl    = -0.5*bsxfun(@rdivide,tmp2,h);
        j0    = dl'*dl/t;
        cov_0 = inv(j0);                                                      

        disp('Analytical results');
        disp( diag(cov_0)'./ diag(cov_qmle)' );

         
        % GARCH(1,1) t-distribution: OPG covariance
        ops   = optimset( 'LargeScale','off','Display','off');  
        start = [a1 b1 nu(k)];
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
% Likelihood function for a GARCH-t(1,1) model
%-------------------------------------------------------------------------
function loglt = lnlt( b,y )
    
    v  = y;                                                                     
    s2 = recserar(1-b(1)-b(2) + b(1)*trimr([0.0; v.^2],0,1),std(y)^2,b(2));     
    nu = abs(b(3));                                             
    z  = v./sqrt(s2);                                                            
    
    const = gamma( (nu+1)/2 ) / ( sqrt(pi*(nu-2)) * gamma( nu/2 ) );           
    loglt = log(const) - 0.5*log(s2) - 0.5*(nu+1)*log( 1 + (z.^2)/(nu-2) );         

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
