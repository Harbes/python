%=========================================================================
%
%   Estimate an LSTAR model of Australian unemployment rate based on the 
%   simple specification to be found in Skalin and Terasvirta ( )
%
%   The delay variable is Delta_4 u
%=========================================================================

function nlm_ozurate( )

    clear all
    clc
    
    % Load data and set up the required lags
    %[data, txt]= xlsread('ausu.xls');
    %load ausu.mat
    load ausu.mat
      
    % Generate LSTAR variables    
    dy      = trimr(data,1,0) - trimr(data,0,1);    % dy(t)
    y_lag   = trimr(data,0,1);                      % y(t-1)
    dy12    = trimr(data,12,0) - trimr(data,0,12);  % y(t)-y(t-12)
    dy12lag = trimr(dy12,0,1);
    
    % Adjust variables to have same length as dy12lag
    dy    = trimr(dy,12,0);
    y_lag = trimr(y_lag,12,0);

    % Create data matrix 
    t = length(dy);     
    
    % Estimate the linear model
    xlin = [ ones( t,1 ) y_lag ];
    bols = xlin\dy;
    lr_linear = -bols(1)/bols(2);
    
     
    % Create elements for the nonlinear model
%     [ r,c ]       = size( xlin );
%     xnln          = (repmat( x(:,end),1,c ).^3).*xlin;
%     [ b,bint,u1 ] = regress( u0,[ xlin xnln ] );
%     RSS1          = sum( u1.^2);
% %     
% %     % Test for nonlinearity
%     test  = t*( (RSS0-RSS1)/RSS1);
%     disp( [test 1-chi2cdf( test,c )] );
%     
    % Optimization
    x       = [ dy y_lag dy12lag ];
    pstart  = 0.1*ones(6,1);
    
    options = optimset( 'Display',             'iter', ...
                        'MaxIter',              2000,   ...
                        'MaxFunEvals',          4000 );   
                    
    phat    = fminunc( @(p) logl(p,x),pstart,options );        
    se      = hess_std( @loglt,phat,x );
    
    disp( [phat se phat./se] );
    
    lr_low = -phat(1)/phat(2);
    lr_high = -(phat(1)+phat(3))/(phat(2)+phat(4));

    disp('Long-run mean unemployment in linear state');
    disp( lr_linear );
    disp('Long-run mean unemployment in low state');
    disp( lr_low );
    disp('Long-run mean unemployment in high state');
    disp( lr_high );

    
            
end

%=========================================================================
%
%   Wrapper function
%
%=========================================================================
function f = logl( p,data )

    f = - sum( loglt( p,data ) );

end


%=========================================================================
%
%   Returns concentrated log-likelihood at each observation
%
%=========================================================================

function f = loglt( p,x )

    y    = x(:,1); 
    ylag = x(:,2);
    st   = x(:,3);
    
    % Transition function
    gam = abs( p(end-1) );
    thr = abs(p(end));
    tmp = (st - thr)/std(st);
    Gt  = 1./( 1+exp( -gam*tmp ) );
 
    % Compute errors
    tmp1 = p(1)+p(2)*ylag; 
    tmp2 = p(3)+p(4)*ylag;
    ut   = y - (tmp1+tmp2.*Gt);
    
    % Concentrate sigma
    s2 = std( ut  );     
   
    % Likelihood function
    f = - 0.5*log(2*pi) - 0.5*log(s2) - 0.5*ut.^2/s2;
    
end

%------------------------------------------------------------------------- 
%   Computes std errors based on Hessian
%------------------------------------------------------------------------- 
function seH = hess_std( llik,bhat,data )

    % Based on the Hessian
    H   = numhess( llik,bhat,data );
    seH = sqrt( diag( inv( H ) ) );    
end

%------------------------------------------------------------------------- 
%   Computes std errors based on Hessian
%------------------------------------------------------------------------- 

function seJ = opg_std( llik,bhat,data )

    % Based on the OPG matrix
    G   = numgrad( llik,bhat,data );
    seJ = sqrt( diag( inv( G'*G ) ) );            
end

%------------------------------------------------------------------------- 
%   Computes numerical gradient at each observation
%-------------------------------------------------------------------------
function G = numgrad( f,x,data )

    f0  = feval( f,x,data );             % n by 1
    n   = length( f0 );
    k   = length( x );
    fdf = zeros( n,k );
 
    % Compute step size 
    dx  = sqrt( eps )*( abs( x ) + eps );
    xh  = x + dx;
    dx  = xh - x;    

    % Compute gradient
    xdx = bsxfun( @plus, diag( dx ), x );
    for i=1:k;
        
        fdf(:,i) = feval( f,xdx(:,i),data );
    end

    G0 = repmat( f0, 1, k );                        % n by k        
    G1 = repmat( dx', n, 1 );
    G  = ( fdf-G0 )./G1;

end


%------------------------------------------------------------------------- 
%   Computes finite difference Hessian
%------------------------------------------------------------------------- 

function H = numhess( f,x,data )
    
    k  = length( x );
    fx = -sum( feval( f,x,data ) );
 
    % Compute the stepsize (h)
    dx  = eps.^( 1/4 )*( abs(x) + eps );
    xh  = x + dx;
    dx  = xh - x;
    ee  = diag( dx ); 
    
    % Compute forward and backward steps
    gplus  = zeros(k,1);
    gminus = zeros(k,1);
    
    for i=1:k
        
        gplus(i)  = -sum( feval( f,x+ee(:,i),data ) );
        gminus(i) = -sum( feval( f,x-ee(:,i),data ) );
    end  
    
    H = bsxfun( @plus, gplus, gminus' ); 
    H = bsxfun( @minus, H, 2*fx )./( dx*dx' );    
    
    
end

