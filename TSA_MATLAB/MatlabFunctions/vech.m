% =======================================================================
%
%   Mimic GAUSS vech( ) function by stacking columns of x
%   on and below the diagonal
%
%========================================================================
function y = vech(x)

    [rows,cols] = size(x);
    y = [];
    for i = 1:rows
        y = [y x(i,1:i)];
    end
    y = y';
end
