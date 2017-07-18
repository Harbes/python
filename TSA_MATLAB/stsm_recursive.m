%=========================================================================
%
%   Recursive structural model estimated in four ways
%
%=========================================================================
function stsm_recursive( )

    clear all
    clc

    % Read the data: quarterly US data from Jan-1959 to Dec-1998
    load simsdata.mat
 
    % Define variables
    r    = ytdata(:,1);        
    lex  = log( ytdata(:,2) );
    lcp  = log( ytdata(:,3) );
    lm   = log( ytdata(:,4) );
    lp   = log( ytdata(:,5) );
    lo   = log( ytdata(:,6) );
    sdum = ytdata(:,7:17);
    
    t = length(ytdata);

    % Construct variables for use in VAR
    % interest rate and the annual percentage growth rates in money, price and output
    yvar = [ r   lm   lp   lo ];    
    tmp  = 100*(trimr(yvar(:,2:4),12,0) - trimr(yvar(:,2:4),0,12));
    yvar = [ trimr(yvar(:,1),12,0) tmp]; 
    lags = [trimr(yvar,1,1)   trimr(yvar,0,2) ];

    % Structural equation approach  
    y  = trimr(yvar(:,1),2,0);
    x  = [ones(length(y),1)   lags ];
    a1 = x\y;
    u1 = y - x*a1;

    y  = trimr(yvar(:,2),2,0);
    x  = [ ones(length(y),1)   trimr(yvar(:,1),2,0)   lags ];
    a2 = x\y;
    u2 = y - x*a2;

    y  = trimr(yvar(:,3),2,0);
    x  = [ ones(length(y),1)   trimr(yvar(:,[1 2]),2,0)   lags ];
    a3 = x\y;
    u3 = y - x*a3;

    y  = trimr(yvar(:,4),2,0);
    x  = [ ones(length(y),1)   trimr(yvar(:,[1 2 3]),2,0)   lags ];
    a4 = x\y;
    u4 = y - x*a4;

    % Structural residuals 	and covariance 
    b01 = [ 1           0        0        0;
           -a2(2)       1        0        0;
           -a3(2)    -a3(3)      1        0;
           -a4(2)    -a4(3)    -a4(4)     1 ];
    u  = [u1 u2 u3 u4 ];  			
    d1 = u'*u/length(u); 		 

    % Reduced form approach	 
    y = trimr(yvar,2,0);
    x = [ones(length(yvar)-2,1)   lags ];

    bar = x\y;
    v   = y - x*bar;
    u1 = v(:,1);
    
    a2 = v(:,1)\v(:,2);
    u2 = v(:,2) - v(:,1)*a2;
    
    a3 = v(:,[1 2])\v(:,3);
    u3 = v(:,3) - v(:,[1 2])*a3;
    
    a4 = v(:,[1 2 3])\v(:,4);
    u4 = v(:,4) - v(:,[1 2 3])*a4;

    % Structural residuals 	and covariance 
    b02 = [ 1           0        0        0;
           -a2(1)       1        0        0;
           -a3(1)    -a3(2)      1        0;
           -a4(1)    -a4(2)    -a4(3)     1 ];
    u  = [u1 u2 u3 u4 ];  			
    d2 = u'*u/length(u); 		 

    % Choleski decomposition approach 
    y = trimr(yvar,2,0);
    x = [ ones(length(y),1)   lags ];
    
    bar    = x\y;
    v      = y - x*bar;
    vc     = v'*v/length(v);
    s      = chol(vc)';
    tmp    = diag(s);
    b0inv  = bsxfun(@rdivide,s,tmp');
    b03    = inv(b0inv);
    d3    = zeros(4);
    d3     = diag(tmp.^2); 

    % MLE approach
    ops    = optimset('LargeScale','off','Display','iter');
    theta0 = 0.1*ones(10,1);
    % Converge in one iteration
    %theta0 = [-b01(2,1) ; -b01(3,1) ; -b01(3,2) ; -b01(4,1) ; -b01(4,2) ; -b01(4,3) ; diag(d1) ]; 

    [ theta,fval ] = fminunc(@(theta) neglog(theta,y,v),theta0,ops);
    lf = -fval;
    disp(['Log-likelihood value     = ',num2str(lf) ]);
    disp(['T x Log-likelihood value = ' num2str(length(v)*lf) ]);

    b04 = [   1                0              0        0 ;
             -theta(1)         1              0        0 ;
             -theta(2)     -theta(3)          1        0 ;
             -theta(4)     -theta(5)      -theta(6)    1 ];

    d4 = zeros(4);
    d4  = diag(theta(7:10));

    disp(' ')
    disp( 'B0: structural equation approach' )	 
    disp( b01 )
    disp( 'B0: reduced form approach' )	 
    disp( b02 )
    disp( 'B0: Choleski decomposition approach' )	 
    disp( b03 )
    disp( 'B0: maximum likelihood approach' )	 
    disp( b04 )


    disp(' ')
    disp( 'D: structural equation approach' )	 
    disp( d1 )
    disp( 'D: reduced form approach' )	 
    disp( d2 )
    disp( 'D: Choleski decomposition approach' )	 
    disp( d3 )
    disp( 'D: maximum likelihood approach' )	 
    disp( d4 )
end

%--------------------------- Functions ----------------------------------
% 
%--------------------------------------------------------------------------
% Log-likelihood function for the VAR
%--------------------------------------------------------------------------
function f = neglog(b,y,v)

    t   = length(v);
	n   = size(y,2);
	lnl = zeros(t,1);

 	b0 =   [   1         0          0      0 ;
         	 -b(1)        1          0      0 ;
	         -b(2)     -b(3)         1      0 ;
             -b(4)     -b(5)      -b(6)    1 ];

    % Structural residual variances 
    d  = eye(n);
    d  = diag(abs(b(7:10)));
	vc = inv(b0)*d*inv(b0)';
     
    for i=1:t
		lnl(i) = -0.5*n*log(2*pi) - 0.5*log(det(vc)) - 0.5*v(i,:)*inv(vc)*v(i,:)';  
    end
    f = -mean( lnl );
end
