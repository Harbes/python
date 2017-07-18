%=========================================================================
%
%   Program to estimate a a capital asset pricing model
%   and estimate the excess return on invested wealth
%
%=========================================================================
function lfac_capm( )

    clear all
    clc

    % Load the data: starting March 1990 and ending March 2010
    % 3 = Microsoft, 
    % 4 = Intel, 
    % 5 = Pfizer, 
    % 6 = Exxon, 
    % 7 = Procter & Gamble, 
    % 8 = Boeing 
    % 9 = AT&T, 
    % 10 = General Electric, 
    % 11 = Chevron, 
    % 12 = Bank of America
    load capm.mat

    sp500    = data(:,1);               % S&P500 index                                                                                                            **/
    interest = data(:,2);               % US Treasury constant maturity, p.a., 3 months                                              **/
    price    = data(:,3:end);           % stocks
                       
    % Compute excess asset excess returns
    tmp = 4*(trimr(log(price),1,0) - trimr(log(price),0,1)); 
    y   = bsxfun(@minus,tmp,trimr(interest,1,0)/100);                                

    % Compute excess market return on sp500
    market_sp500 = 4*(trimr(log(sp500),1,0) - trimr(log(sp500),0,1)) - trimr(interest,1,0)/100;                                
    
    t = length(y);
    
    % Estimate parameters with restriction imposed
    flag  = 1;
    start = [ones(20,1) ; 1 ; ones(10,1) ;];
    ops   = optimset('LargeScale','off','Display','iter', ...
                    'MaxFunEvals',Inf,'MaxIter',Inf);
    
    bhat = fminunc(@(b) neglog( b,y,flag ),start,ops );
    
    % Restimate without restrictions to get standard errors
    flag  = 0;
    start = [bhat(1:20) ; tanh(bhat(21)) ; bhat(22:end)] ;

        
   [ bhat,lf,~,~,~,hess ] = fminunc(@(b) neglog( b,y,flag ),start,ops );
   
    disp( ['Log-likelihood function = ',num2str(-lf) ] );
    disp( ' ' )

    % Estimate excess returns on wealth - mean  
    [~,s11] = kfsmooth(bhat,y);
    
    % CAPM regression estimates  
    beta = [ones(t,1) market_sp500]\y;

    disp(' Estimates of Beta ')
    disp('    OLS       KF        KF(scaled) ')
    disp( [beta(2,:)' bhat(11:20) bhat(11:20)*std(s11)/std(market_sp500) ])



end

%
% ------------------------ Functions ------------------------------------%
%
%-------------------------------------------------------------------------
%  Log-likelihood function from Kalman filter
%-------------------------------------------------------------------------
function lf = neglog( b,y,flag )

    % Define matrices for filter
    alpha = b(1:10);      
    Lam   = b(11:20);
    if flag
        Phi = tanh(b(21));
    else
        Phi = b(21);
    end
    R     = diag(b(22:31).^2);
    Q     = diag(length(Phi));     
    
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
    alpha = b(1:10);      
    Lam   = b(11:20);
    if flag
        Phi = tanh(b(21));
    else
        Phi = b(21);
    end
    R     = diag(b(22:31).^2);
    Q     = diag(length(Phi));     
    
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
   


