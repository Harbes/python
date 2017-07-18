%========================================================================
%
%   Estimate a Structural VAR based on long-run restrictions
%
%========================================================================
function svar_longrun( )

    clear all;
    clc;

    % U.S.zero coupon yields Dec 1946 - Feb 2002 (0,1,3,6,9 mth and 10 yrs)
    load mcnew.mat
    yt   = trimr(rt(:,1:3),1,0)-trimr(rt(:,1:3),0,1);

    % Estimate a VAR with one lag
    y  = trimr(yt,1,0);
    x  = [ones(length(y),1)  trimr(yt,0,1) ];
    a  = x\y;
    v  = y - x*a;                            
    vc = v'*v/length(v);
    
    % Lag 1 autoregressive parameters   
    a1 = trimr(a,1,0)';     

    % Esimate SVAR
    opt      = optimset('LargeScale','off','Display','iter', ...
                        'MaxIter',Inf,'MaxFunEvals',Inf);
    theta0   = ones(5,1);
    [b,fval] = fminunc(@(b) neglog(b,v,a1),theta0,opt);
    
    disp(['Log-likelihood value = ', num2str(-fval) ]);
    disp(' ')
    disp('VAR parameter estimates')
    disp( b )
    disp(' ')
    disp('VAR variance-covariance matrix (unconstrained)')
    disp( vc )
    
    % Compute SVAR matrices
    [f,s,vc0] = func(b,a1);
    
    disp('VAR variance-covariance matrix (constrained)')
    disp( vc0 )
    disp(['Determinant of omegav   = ', num2str(det(vc0)) ]);
    disp('I-A1')
    disp( eye(3) - a1  )
    disp('f') 		
    disp( f );
    disp('s') 		
    disp( s );
    
end
%
%-------------------------Functions------------------------------------
%
%-----------------------------------------------------------------------
%      Negative log-likelihood function   
%-----------------------------------------------------------------------
function lf = neglog(b,v,a1)

    [ t,n ] = size(v);
	 
 	f =   [    b(1)     0      0   ;
         	   b(4)    b(2)    0   ;
	           b(5)     0      b(3) ];
  
    s  = (eye(3) - a1)*f;
	vc = s*s';

    lnl=zeros(t,1);
    for i=1:t

		lnl(i) = -0.5*n*log(2*pi) - 0.5*log(det(vc)) ...
                 - 0.5*v(i,:)*inv(vc)*v(i,:)';    
             
     end
     lf = -mean( lnl );
end
%-----------------------------------------------------------------------
%   Return SVAR matrices   
%-----------------------------------------------------------------------
function [f,s,vc0] = func(b,a1)

	 
 	f  =   [    b(1)     0      0   ;
         	    b(4)    b(2)    0   ;
	            b(5)     0      b(3) ]; 
    s   = (eye(3) - a1)*f;
	vc0 = s*s';


end


