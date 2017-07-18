%=========================================================================
%
%    Kalman Filter implementation of the Hodrick Prescott Filter 
%
%=========================================================================

function lfac_hp( )
    
    clear all
    clc
    
    % Load quarterly US data for the period 1940:1 to 2000:4 (T = 244)
    load lfac_usgdp.mat    
    %[RGDP, txt]= xlsread('usgdp.xlsx');
    
    y = log(RGDP)-mean(log(RGDP));
    y = y*100;
    t = length(y);
    
    % Estimate unconstrained model
    start = [ 1 ; 1] ;
    ops   = optimset('LargeScale','off','Display','off');
          
    [bhat,lf1] = fminunc(@(b) neglog(b,y),start,ops );
    
    lf1 = -lf1;
    
    disp('Parameter Estimates')
    disp( bhat )
    
    % Compute likelihood from the Kalman filter with HP restrictions
    [ fac,lf0 ] = hpsmooth( y );
    
    trend_kf = fac(:,1);
    
    lf0 = -lf0;
    
    disp(['Log-likelihood function (unconstrained)   = ',num2str(lf1) ]); 
    disp(['Log-likelihood function (constrained)     = ',num2str(lf0) ]); 
    
    lr = -2*t*(lf0 - lf1);

    disp(['LR statistic    = ',num2str(lr) ]);
    disp(['p-value         = ',num2str(1-chi2cdf(lr,1)) ]);



    % Get trend component from HP filter and plot results
    
    trend_hp = HPfilter( y,1600 );

    %trend_hp = hpfilter( y,1600 ); % Matlab version gives same results!


    
    %**********************************************************************
    %***
    %***     Generate graph
    %***
    %**********************************************************************
    %dt = datenum(txt,'yyyy-q');
    dt = 1:54;
    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;

    %plot( [log(RGDP),Pc2 ])
    plot( dt(1:54),y(1:54),'--k',dt(1:54),trend_hp(1:54),'-k','LineWidth',1 );
    datetick('x','yyyy')
    ylabel('Log real US GDP $\times 100$');
    xlabel('t');
    box off
    axis tight

       
end

%--------------------------- Functions ----------------------------------
% 
%-------------------------------------------------------------------------
%       Univariate Hodrick Prescott Filter in Kalman Filter Form
%       Unrestricted so parameters can be estimated.
%-------------------------------------------------------------------------
function lf = neglog(b,y)

    % Set up filter 
    
    Lam  = [1 0];
    Phi  = [ 1 1; 0 1];

    R = b(1)^2;
    Q = [ 0 0; 0 b(2)^2];

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
    pt = eye(k)*1000;   

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
%-------------------------------------------------------------------------
%       Univariate Hodrick Prescott Filter in Kalman Filter Form
%-------------------------------------------------------------------------
function [lfac,lf] = hpsmooth( y )

    % Set up filter
    sigu = 1;
    sigz = sqrt(1/1600);
    
    Lam = [1 0];
    Phi = [ 1 1; 0 1];

    % R
    R = sigu;
    
    % Q
    Q = [ 0 0; 0 sigz];
	
    % Allocate arrays
    [ t,n ] = size(y);
    lnl     = zeros(t,1);
    k       = size(Q,1);           
    s10     = zeros( t,k );    % st|t-1
    s11     = zeros( t,k );    % st|t	
    p10     = zeros( k,k,t );
    p11     = zeros( k,k,t );
    ss      = zeros(t,k);
  	   
    % Initialisation following Harvey and Hamilton  
    st = zeros(k,1);
    pt = eye(k)*1000;  
    s10(1,:) = st; 
    p10(:,:,1) = pt;
 
    mt = Lam*st;
    vt = Lam*pt*Lam' + R;
    ut = y(1,:)' - mt;
    lnl(1) = - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5*ut'*inv(vt)*ut;
    
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
                       
        lnl(i) = - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5*ut'*inv(vt)*ut;
        
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

   lfac = ss;   
   lf   = -mean(lnl);

end
%-------------------------------------------------------------------------
%  Univariate Hodrick Prescott Filter in Original Form
%-------------------------------------------------------------------------
function HPtrend = HPfilter( y,lambda )

    [t,m] = size (y);

    if t < m
        y = y';     
        t = m;
    end
    
    % Setting up
    td   = t - 4;
    tmp  = zeros(t,1);
    tmpd = zeros(td,1);

    % Set up the first and two rows of the filter matrix F
    r1 = tmp;
    r1(1) = 1.0;  r1(2) = -2.0; r1(3) = 1.0;
    
    r2 = tmp;
    r2(1) = -2.0; r2(2) = 5.0;  r2(3) = -4.0; r2(4) = 1.0;
    
    
    % Create the diagonals for the filter matrix F
    tmp  = ones(t,1);
    tmpd = ones(td,1);

    D = [tmpd -4.0*tmpd 6.0*tmpd -4.0*tmpd tmpd];
    
    % Construct the filter matrix F
    F = [ r1 r2 spdiags( D,-4:0,t,td ) flipud(r2) flipud(r1) ];
    
    
    HPtrend =inv(lambda*F+eye(t))*y;
end

