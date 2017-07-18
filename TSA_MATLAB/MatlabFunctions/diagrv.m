%========================================================================
%
%   Puts given vector (v) onto the main diagonal of a square matrix (x)
%   Reproduces the GAUSS command diagrv
%
%========================================================================
function y = diagrv(x,v)

y = x - diag(diag(x)) + diag(v);

end

