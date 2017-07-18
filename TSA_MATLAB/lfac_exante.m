%=========================================================================
%
%   Program to to estimate ex ante real interest rates  
%   from ex post real interest rates
%
%=========================================================================
function lfac_exante( )

    clear all
    clc


    % Load data starting Jan 1971 and ending December 2009
    load exante.mat
    
    % Compute ex post real interest rate    
    inflation = 1200*(trimr(log(price),1,0) - trimr(log(price),0,1));           
    y = trimr(interest,1,0) - inflation;                                

    t = length(y);
    
    % Estimate parameters with restriction imposed
    flag  = 1;
    start = [mean(y) ; 0.5 ; 1 ; std(y) ];
    ops   = optimset('LargeScale','off','Display','off');
    
    bhat = fminunc(@(b) neglog( b,y,flag ),start,ops );
    
    % Restimate without restrictions to get standard errors
    flag  = 0;
    start = [bhat(1); tanh(bhat(2)) ; bhat(3:4) ];
        
    [ bhat,lf,~,~,~,hess ] = fminunc(@(b) neglog( b,y,flag ),start,ops );

    disp( ['Log-likelihood function = ',num2str(lf) ] );
    disp( ' ' )
    disp( [' Mean of ex post real interest rate mean     = ',num2str(mean(y)) ]);
    disp( [' Variance of ex post real interest rate mean = ',num2str(std(y)^2) ]);

    disp( ' ' )   
    disp( [' Mean of ex ante real interest rate mean     = ',num2str(bhat(1)) ]);
    disp( [' Variance of ex ante real interest rate mean = ',num2str(bhat(4)^2/(1 - bhat(2)^2)) ]);
    
    % Real ex ante interest rate
    [~,s_update] = kfsmooth(bhat,y);
    
    % Plot leaves off addition of bhat(1) so you can see two series
    exante = s_update;% + bhat(1); 
    figure(1)
    plot(seqa(1971+2/12,1/12,t),[exante y]);
    title('Ex Ante Interest Rate')
        
end
%
% ------------------------ Functions ------------------------------------%
%
%-------------------------------------------------------------------------
%  Log-likelihood function from Kalman filter
%-------------------------------------------------------------------------
function lf = neglog( b,y,flag )

    % Define matrices for filter
    alpha = b(1);      
    Lam   = 1;
    if flag
        Phi = tanh(b(2));
    else
        Phi = b(2);
    end
    R     = b(3)^2;
    Q     = b(4)^2;     
    
    lf    = -mean( lnlt(y,Phi,Lam,R,Q,alpha) ); 
	
end
%--------------------------------------------------------------------------
% Multivariate Kalman filter
%--------------------------------------------------------------------------
function lnl= lnlt(y,Phi,Lam,R,Q,alpha)


    % Allocate arrays
    [ t,n ]   = size(y);
    k         = size(Q,1);
    lnl       = zeros(t,1);
    
	% Recursions of the Kalman Filter
    % Initialisation following Harvey and Hamilton  
    st = zeros(k,1);
    pt = reshape(inv(eye(k^2) - kron(Phi,Phi))*Q(:),k,k)';   

    mt = Lam*st + alpha;
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
        mt = Lam*st + alpha;
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
%--------------------------------------------------------------------------
%   Extract smoothed factor
%--------------------------------------------------------------------------
function [fac,s11] = kfsmooth(b,y)

    % Unpack the parameter vector
    alpha = b(1);      
    Lam   = 1;
    if flag
        Phi = tanh(b(2));
    else
        Phi = b(2);
    end
    R     = b(3)^2;
    Q     = b(4)^2;     
    
    % Allocate arrays
    [ t,n ] = size(y);
    k       = size(Q,1);           
    s10     = zeros( t,k );    % st|t-1
    s11     = zeros( t,k );    % st|t	
    p10     = zeros( k,k,t );
    p11     = zeros( k,k,t );
    ss      = zeros(t,1);
  	   
    % Initialisation following Harvey and Hamilton  
    st = zeros(k,1);
    pt = eye(k)*0.1;  
    s10(1,:) = st; 
    p10(:,:,1) = pt;
 
    mt = Lam*st + alpha;
    vt = Lam*pt*Lam' + R;
    ut = y(1,:)' - mt;
    
    Kal = pt*Lam'*inv(vt);
    s0 = st + Kal*ut;
    p0 = pt - Kal*Lam*pt;
    s11(1,:) = s0;
    p11(:,:,1) = p0;


    % Main loop over observations

    for i = 2:t
        
	    % Prediction 
        st = Phi*s0;                     
        pt = Phi*p0*Phi' + Q;	    
        s10(i,:) = st; 
        p10(:,:,i) = pt;
      
        % Observation
        mt = Lam*st + alpha;
        vt = Lam*pt*Lam' + R;
        ut = y(i,:)' - mt;
                       
		% Update    			
        Kal = pt*Lam'*inv(vt);
        s0 = st + Kal*ut;
        p0 = pt - Kal*Lam*pt;
        s11(i,:) = s0;
        p11(:,:,i) = p0;
        
    end
    
    % Now smooth the factor    
    ss = s11;
    for j = 1:t-1
       
        Jt      = p11(:,:,t-j)*Phi'*inv(p10(:,:,t-j));     
        ss(t-j,:) = s11(t-j,:)' + Jt*( ss(t-j+1,:)' - s10(t-j+1,:)' ); 
    
    end
    fac = ss;   

end
   

