%==========================================================================
%
%   	Program to estimate a symmetric BEKK MGARCH model of US yields
%       with conditionsl multivariate Student t disturbance
%
%==========================================================================
function mgarch_student( )

    clear all
    clc

    % Load data
    load yields_us

    % Choose variables
    r  = rdata(:,[ 1 2 ] );                                                  
    y  = 100*(trimr( r,1,0 ) - trimr( r,0,1 )); 
    t  = length( y );
 
    % Estimate the BEKK Aymmetric MGARCH(1,1) model             
    start = [   0.8062958419864796 
                0.2816730578640683 
                0.4297482870507287 
                0.2685150842768220 
               -0.0075466373201676 
                0.1615508008596290 
                0.9269147664364472 
                0.0083884990798515 
                0.9724763092352805 
               -0.1336772926059021 
               -0.0906508454299734 
                8.3673242655198212 ];

    ops = optimset( 'LargeScale','off','Display','iter' );
 
    [theta,~ ] = fminunc( @(b) neglog(b,y),start,ops); 
    
    disp('BEKK parameters')
    disp(theta)
end
    
%
%--------------------------- Functions ----------------------------------
% 
%--------------------------------------------------------------------------
% Log-liklihood function for symmetric BEKK model
%--------------------------------------------------------------------------
function lf = neglog( b,y )

    [ t,n ] = size( y );
    f      = zeros( t,1 );

    h = cov(y);     % Initialise  conditional covariance matrix
    v = std(y);     % Initialise the disturbance vector
    v = v';
    
    c  = [ b(1)    0.0  ;     
           b(2)    b(3) ];

    a =  [ b(4)    b(5) ; 
           b(5)    b(6) ];

    d =  [ b(7)    b(8) ;
           b(8)    b(9) ];
       
    nu = b(12);

    for i = 1:t
        
        m    = b([10 11]);                   % Update conditional mean
        v    = y(i,:)'- m;                   % Update residuals                                

       const = gamma( (nu+1)/2 ) / ( sqrt(pi*(nu-2)) * gamma( nu/2 ) );                    
       f(i)  = log(const) - 0.5*log(det(h)) - 0.5*(nu+1)*log( 1 + v'*inv(h)*v/(nu-2) );     
              
        h    = c*c' + a*(v*v')*a' + d*h*d';  % Update conditional covariance matrix 


    end
    lf = -mean( f );
    
end

