%------------------------------------------------------------------------- 
%   Computes finite difference Hessian
%------------------------------------------------------------------------- 

function H = numhess( f,x,varargin )
    
    k  = length( x );
    f0 = feval( f, x, varargin{:} );
 
    % Compute the stepsize (h)
    dx  = eps.^( 1/3 )*( abs(x) + eps );
    xh  = x + dx;
    dx  = xh - x;
    ee  = diag( dx ); 
 
    % Compute forward and backward steps
    fplus  = zeros( k,1 );
    
    for i=1:k
        
        fplus(i)  = feval( f, x+ee(:,i), varargin{:} );
    end  
    
    H = zeros( k );
    for j = 1:k
       
        for l = 1:j;
            
            H(j,l) = feval( f, x+ee(:,j)+ee(:,l), varargin{:} ); 
        end
    end
    H = H + tril( H,-1 )';    
    
    fpp = bsxfun( @plus, fplus, fplus' ); 
    H   = H - fpp + f0;
    H   = H./( dx*dx' );    
    
end

