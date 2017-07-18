%=========================================================================
%
%   Program to estimate a one factor model of the term structure
%   Alternative formulation of the Kalman filter
%
%=========================================================================
function lfac_term_adj( )

    clear all;
    clc;

    % Read the data:
    % US daily zero coupon yields starting 10/04/1988 and ending 12/28/2001
    % 	The variables are:
    %
    %		1.	tcm3m
    %		2.	tcm1y
    %		3.	tcm3y
    %		4.	tcm5y
    %		5.	tcm7y
    %		6.	tcm10y

    load lfac_usdata.mat
    interest = usdata;

    %	Choose the variables				
    r = interest(:,[6 5 4 3 2 1]);

    %	Rescale the variables (used for numerical precision when using mle)		
    y = r*100;                   
    y = bsxfun(@minus,y,mean(y));
    t = length(y);

        % Estimate parameters with restriction imposed
    flag  = 1;   
    start = [   5.5359 
                5.8476   
                6.2928    
                7.0033   
                7.6360   
                7.3535   
               50.4656   
               36.6785  
               26.5888   
                6.2867   
               43.6677  
               65.7365  
                0.9993];
    ops   = optimset('LargeScale','off','Display','iter', ...
                      'MaxFunEvals',Inf,'MaxIter',Inf);
    
    bhat = fminunc(@(b) neglog( b,y,flag ),start,ops );
    
    % Restimate without restrictions to get standard errors
    flag  = 0;
    start = [bhat(1:12); tanh(bhat(13)) ];
        
    [ bhat,lf,aa,aa,aa,hess ] = fminunc(@(b) neglog( b,y,flag ),start,ops );

    vc = (1/t)*inv(hess);

    disp( ['Log-likelihood function = ',num2str(-lf) ] );
    disp( ' ' )
    disp('  Params   Std. Errors   ')
    disp( [bhat  sqrt(diag(vc)) ] );


end
%
% ------------------------ Functions ------------------------------------%
%
%-------------------------------------------------------------------------
%  Log-likelihood function from Kalman filter
%-------------------------------------------------------------------------
function lf = neglog( b,y,flag )

    % Define matrices for filter
    Lam = [b(1:6)   diag(abs(b(7:12))) ];
 
    if flag
        Phi = diag( [tanh(b(13)) ; zeros(6,1)] );
    else
        Phi = diag( [b(13) ; zeros(6,1)] );
    end
    R     = zeros(6,6);
    Q     = eye(7);     
    
    lf    = -mean( lnlt(y,Phi,Lam,R,Q) ); 
	
end
%--------------------------------------------------------------------------
% Multivariate Kalman filter
%--------------------------------------------------------------------------
function lnl= lnlt(y,Phi,Lam,R,Q)


    % Allocate arrays
    [ t,n ]   = size(y);
    k         = size(Q,1);
    lnl       = zeros(t,1);
    
	% Recursions of the Kalman Filter
    % Initialisation following Harvey and Hamilton  
    st = zeros(k,1);
    pt = eye(k);
    pt(1,1) = 0.1;    % Impose prior on first factor
    
    mt = Lam*st;
    vt = Lam*pt*Lam' + R;
    ut = y(1,:)' - mt;

    lnl(1) = - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5*ut'*inv(vt)*ut;
    
    Kal = pt*Lam'*inv(vt);
    s0 = st + Kal*ut;
    p0 = pt - Kal*Lam*pt;

    % Main loop over observations

    for i = 2:t
        
	    % Prediction 
        st = Phi*s0;                     
        pt = Phi*p0*Phi' + Q;	        
              
        % Observation
        mt = Lam*st;
        vt = Lam*pt*Lam' + R;
        ut = y(i,:)' - mt;
       
        % Construct log-likelihood function
        lnl(i) = - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5*ut'*inv(vt)*ut;
                
		% Update    			
        Kal = pt*Lam'*inv(vt);
        s0 = st + Kal*ut;
        p0 = pt - Kal*Lam*pt;
        
    end
end
