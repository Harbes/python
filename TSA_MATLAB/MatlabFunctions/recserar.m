%=========================================================================
% Computes a vector of autoregressive recursive series
%
%   Inputs: 
%             x  = a matrix of dimensions (n,k) (exogenous variables)
%             y0 = a matrix of dimensions (p,k) (starting values)
%             a  = a matrix of dimensions (p,k) (lag parameters)
%
%  Mimics the Gauss routine. 
%=========================================================================



function y = recserar(x,y0,a)

     
    [rx cx] = size(x);
    [ry cy] = size(y0);
    [ra ca] = size(a);

    if (nargin ~= 3) || (cx ~= cy) || (cx ~= ca) || (ry ~= ra) 
        error('Check function inputs');
    end
 
    y=zeros(rx,cx);

    for j=1:ry;
        y(j,:) = y0(j,:);
    end;
    for j=(ry+1):rx;
        y(j,:)=x(j,:);
            for k=1:ry;
                y(j,:)=y(j,:)+a(k,:).*y(j-k,:);
            end;
     end;
end
