%==========================================================================
%
%   Program to estimate a international capital asset pricing model 
%   with time-varying beta risk using the symmetric BEKK specification.
%
%==========================================================================
function mgarch_icapm( )

    clear all
    clc

    % Daily data starting 3-Feb-1988 to 29-Dec-1995
	%	1.	date
	%	2.	interest rate (risk free)
	%	3.	nyse returns (NYSE stock index)
	%	4.	world returns (from MSCI index)

    load icapm.mat
    
    r = 100*(rdata(:,3) - rdata(:,2));      % NYSE excess return          
    m = 100*(rdata(:,4) - rdata(:,2));      % World excess return         
    t = length(r);                                    
    y = [ r  m ];
    
    
    % Constant beta risk estimate (based on OLS)      
    x     = [ ones(t,1)   m ];
    b_ols = x\r;

    disp(['Constant estimate of beta risk = ', num2str(b_ols(2)) ]);

    % Estimate the BEKK Symmetric MGARCH(1,1) model             
    start = [   0.065125621828620
               -0.082834944191625
                0.054312771634006
                0.141289906762252
               -0.037479056981814
                0.257666483252922
                0.974472477018916
                0.024314746700112
                0.946834227119610
                0.063159015564876
                0.040511737867333 ];
    ops = optimset( 'LargeScale','off','Display','iter');
 
    [theta,~ ] = fminunc( @(b) neglog(b,y),start,ops); 
         
    % Compute beta risk at optimal parameters
    [ ~,br ] = neglog( theta,y );
    
    % Plot conditional beta risk
    figure(1)
    plot( seqa(1988+2/12,1/257,t) , br );
    title('Conditional Beta Risk');
    xlabel('t');
    ylabel('$\beta_t$')
    
end
%
%--------------------------- Functions ----------------------------------
% 
%--------------------------------------------------------------------------
% Log-likelihood function for symmetric BEKK model
%--------------------------------------------------------------------------
function [lf,br] = neglog( b,y )

    [ t,n ] = size( y );
    f      = zeros( t,1 );
    br     = zeros( t,1 );

    h = cov(y);     % Initialise  conditional covariance matrix
    u = std(y);     % Initialise the disturbance vector
    u = u';
    
    c  = [ b(1)    0.0  ;     
           b(2)    b(3) ];

    a =  [ b(4)    b(5) ; 
           b(5)    b(6) ];

    d =  [ b(7)    b(8) ;
           b(8)    b(9) ];

    for i = 1:t
               
        h    = c*c' + a*(u*u')*a' + d*h*d';  % Update conditional covariance matrix 

        m    = b([10 11]);                   % Update conditional mean
        u    = y(i,:)'- m;                   % Update residuals                                
        f(i) = -0.5*n*log(2*pi) - 0.5*log(det(h)) - 0.5*u'*inv(h)*u;       

        br(i) = h(1,2)/h(2,2);

    end
    lf = -mean( f );
    
end





