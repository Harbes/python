%------------------------------------------------------------------------- 
%   Computes numerical gradient at each observation
%-------------------------------------------------------------------------
function G = numgrad( f,x,varargin )

    f0  = feval( f,x,varargin{:} );             % n by 1
    n   = length( f0 );
    k   = length( x );
    fdf = zeros( n,k );
 
    % Compute step size 
    dx      = sqrt( eps )*( abs( x ) + eps );
    xh      = x + dx;
    dx      = xh - x;    
    ind     = dx < sqrt(eps);
    dx(ind) = sqrt(eps);

    % Compute gradient
    xdx = bsxfun( @plus, diag( dx ), x );
    for i=1:k;
        
        fdf(:,i) = feval( f, xdx(:,i), varargin{:} );
    end

    G0 = repmat( f0, 1, k );                        % n by k        
    G1 = repmat( dx', n, 1 );
    G  = ( fdf-G0 )./G1;

end
