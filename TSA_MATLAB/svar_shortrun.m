%========================================================================
%
%   Estimate a Structural VAR based on short-run restrictions
%
%========================================================================
function svar_shortrun( )

    clear all;
    clc;

    % U.S.zero coupon yields Dec 1946 - Feb 2002 (0,1,3,6,9 mth and 10 yrs)
    load mcnew.mat
    yt   = rt(:,1:3);

    % Estimate a VAR with one lag
    y  = trimr(yt,1,0);
    x  = [ones(length(y),1)  trimr(yt,0,1) ];
    a  = x\y;
    v  = y - x*a;                            
    vc = v'*v/length(v);

    % Esimate SVAR
    opt      = optimset('LargeScale','off','Display','iter');
    theta0   = 0.1*ones(5,1);
    [b,fval] = fminunc(@(b) neglog(b,v),theta0,opt);
    
    disp(['Log-likelihood value = ', num2str(-fval) ]);
    disp(' ')
    disp('VAR parameter estimates')
    disp( b )
    disp(' ')
    disp('VAR variance-covariance matrix (unconstrained)')
    disp( vc )
    
    % Compute SVAR matrices
    [b0,d,s,vc0] = func(b,v);
    
    disp('VAR variance-covariance matrix (constrained)')
    disp( vc0 )
    disp(['Determinant of omegav   = ', num2str(det(vc0)) ]);
    disp('b0')
    disp( b0 )
    disp('d') 		
    disp( d );
    disp('s') 		
    disp( s );
    
end
%
%-------------------------Functions------------------------------------
%
%-----------------------------------------------------------------------
%      Negative log-likelihood function   
%-----------------------------------------------------------------------
function lf = neglog(b,v)

    [ t,n ] = size(v);
	 
 	b0 =   [   1       0      0 ;
         	  -b(1)    1      0 ;
	          -b(2)    0      1 ];

    % Structural residual variances 
	d  = diagrv(eye(n),abs(b(3:5)));		           
    s  = inv(b0)*sqrt(d);
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
function [b0,d,s,vc0] = func(b,v)

    [ ~,n ] = size(v);
	 
 	b0 =   [   1       0      0 ;
         	  -b(1)    1      0 ;
	          -b(2)    0      1 ];

    % Structural residual variances 
	d   = diagrv(eye(n),abs(b(3:5)));		           
    s   = inv(b0)*sqrt(d);
	vc0 = s*s';

end


