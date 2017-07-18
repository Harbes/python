% =========================================================================
%
% Reshape a matrix to agree with GAUSS reshape command
%
% =========================================================================

function X = reshapeg(Y,rows,cols)

         tmp = reshape(Y',cols,rows);
         X   = tmp';

end