%=========================================================================
%
% Column stack a matrix for agree with GAUSS vec command
%
%=========================================================================
function y = vec(x)

    n = numel(x) ;
    y = reshape( x,[n,1] );

end