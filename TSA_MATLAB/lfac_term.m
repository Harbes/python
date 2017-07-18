%========================================================================
%
%		Estimate a one factor model of the term structure
%
%========================================================================
function lfac_term( )

    clear all;
    clc;

    % Read the data, rearrange and scale data
    % US daily zero coupon yields starting 10/04/1988 and ending 12/28/2001
    %		1.	tcm10y
    %		2.	tcm7y
    %		3.	tcm5y
    %		4.	tcm3y
    %		5.	tcm1y
    %		6.	tcm3m
    
    load lfac_usdata.mat;
    rt = usdata*100;
    yt = bsxfun(@minus,rt,mean(rt)); 
    t  = length(yt);

    % Estimate the model by MLE	
    start = [  7.353852219681254 
                7.636371949157256 
                7.003563287217877 
                6.293071514715054 
                5.847811722064529 
                5.536175149488106 
                65.74336399439356 
                43.66727769088171 
                6.286986479990421 
                26.59070018994360 
                36.68045957160562 
                50.46637407916106 
                3.979446498229009 ];
 
            
    % Estimate model
    opt   = optimset('LargeScale','off','Display','iter', ...
                    'MaxFunEvals',Inf,'MaxIter',Inf);
    
    [theta1,fval1] = fminunc(@(b) neglog(b,yt,1),start,opt);
    
    fval1 = -fval1;
    theta1(end) = tanh(theta1(end));
    
   
    disp( ' ' )
    disp('  Params   ')
    disp( theta1   );

    
    % Extract and plot the smoothed factor
    fac = kfsmooth(theta1,yt);
    figure(1)
    plot(1:length(fac),fac);
    
    % Estimate restricted model
    start = [
                    6.862810483385435 
                    66.30255435166845 
                    45.85160110029478 
                    6.589040449849120 
                    29.53625724264360 
                    42.59633524199530 
                    57.62141009765852 
                    3.988334426025099 ];
    [theta0,fval0] = fminunc(@(b) neglogc(b,yt,1),start,opt);
    fval0 = -fval0;
    % Likelihood ratio test of restrictions
    lr  = -2*t*(fval0 - fval1);
    dof = length(theta1) - length(theta0);

    disp(['Log-likelihood (unrestricted) = ',num2str(fval1) ]);
    disp(['Log-likelihood (restricted)   = ',num2str(fval0) ]);
    disp(['LR statistic                  = ',num2str(lr) ]);
    disp(['p-value                       = ',num2str(1-chi2cdf(lr,dof)) ]);


end
%--------------------------- Functions ----------------------------------
% 
%--------------------------------------------------------------------------
%   Wrapper function to set up and call the Kalman filter
%--------------------------------------------------------------------------
function lf = neglog(b,y,flag)

    % Unpack the parameter vector
    Lam = b(1:6);
    if flag
        Phi = tanh(b(13));
    else
        Phi = b(13);
    end

    %R = eye(6);
    R = diag(b(7:12).^2);
    Q = 1;
    
    lf = -mean( kalman(y,Phi,Lam,R,Q) );
   
end
%--------------------------------------------------------------------------
%   Wrapper function to set up and call the Kalman filter
%--------------------------------------------------------------------------
function lf = neglogc(b,y,flag)

    % Unpack the parameter vector
    Lam = [ b(1) ;
            b(1) ;
            b(1) ;
            b(1) ;
            b(1) ;
            b(1) ];
    
        if flag
        Phi = tanh(b(8));
    else
        Phi = b(8);
    end

    %R = eye(6);
    R = diag(b(2:7).^2);
    Q = 1;
    
    lf = -mean( kalman(y,Phi,Lam,R,Q) );
   
end

%--------------------------------------------------------------------------
% Kalman filter
%--------------------------------------------------------------------------
function lnl = kalman(y,Phi,Lam,R,Q)


    % Allocate arrays
    [ t,n ]   = size(y);
    k         = size(Q,1);
    lnl       = zeros(t,1);
  	
    % Initialisation following Harvey and Hamilton  
    st = zeros(k,1);
    pt = eye(k)*0.1;   

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

%--------------------------------------------------------------------------
%   Extract smoothed factor
%--------------------------------------------------------------------------
function fac = kfsmooth(b,y)

    % Unpack the parameter vector
    Lam = b(1:6);
    Phi = b(13);
    R   = diag(b(7:12).^2);
    Q   = 1;
    
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
 
    mt = Lam*st;
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
        mt = Lam*st;
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
   

