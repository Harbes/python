%=========================================================================
%
%   Kalman filter program to estimate a business cycle model for Australia
%
%=========================================================================
function lfac_bcycle( )

    clear all
    clc

    % Load the data: 
    % Australian monthly business cycle data Jan 1980 - September 2009
    load lfac_bcycle.mat


    coincident   = data(:,1);        % Coincident index                                      
    gdp          = data(:,2);        % GDP (interpolated)               
    unemployment = data(:,3);        % Unemployment rate (p.a.)                             
    employment   = data(:,4);        % Employment                                           
    sales        = data(:,5);        % Retail sales                                          
    income       = data(:,6);        % Household income (interpolated)
    production   = data(:,7);        % Industrial production                                

    ind = [gdp unemployment  employment sales income production];

    % Transformed indicators 
    y = [ 100*(trimr(log(ind(:,1)),12,0) - trimr(log(ind(:,1)),0,12)) ...
               trimr(ind(:,2),12,0) - trimr(ind(:,2),0,12)            ...
          100*(trimr(log(ind(:,3:6)),12,0) - trimr(log(ind(:,3:6)),0,12)) ];                    
    y = bsxfun(@minus,y,mean(y));
    t = length(y);


    % Estimate parameters 
    %start = [ones(6,1) ; 0.9 ; -0.5 ; ones(6,1) ];   
    start  = [   0.0917      
                -0.0584        
                 0.0836         
                 0.0644         
                 0.0227         
                 0.0934         
                 1.8952         
                -0.9125         
                 0.9376         
                -0.2106        
                 0.7740         
                 2.8910         
                 3.0981         
                 3.8242   ];      

    ops   = optimset('LargeScale','off','Display','iter', ...
                    'MaxFunEvals',Inf,'MaxIter',Inf);
    
    [bhat,lf,~,~,~,hess] = fminunc(@(b) neglog( b,y ),start,ops );
    
    lf = -lf;
    vc = (1/t)*inv(hess);
    disp(['Log-likelihood function    = ',num2str(lf) ])
    disp( '  Params    Std Errors ')
    disp( [bhat sqrt(diag(vc))] )

    % Compute business cycle based on coincident index  
    bc = 100*(trimr(log(coincident),12,0) - trimr(log(coincident),0,12));      
    bc = bc - mean(bc);

    % Compute and scale the smoothed factor estimate 
    fac = kfsmooth(bhat,y);
    bm  = fac(:,1);
    bm  = bm*std(bc)/std(bm);

    figure(1)
    plot(seqa(1981+1/12,1/12,t),[bm bc]);
    
end
%--------------------------- Functions ----------------------------------
% 
%--------------------------------------------------------------------------
%   Wrapper function to set up and call the Kalman filter
%--------------------------------------------------------------------------
function lf = neglog(b,y)

    Lam = [ b(1:6) zeros(6,1)];
    Phi = [ b(7)  b(8) ;
            1     0   ];
    
    R   = diag(b(9:14).^2);
    Q = diag([1;zeros(length(Phi)-1,1)]);
    
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
    Lam = [ b(1:6) zeros(6,1)];
    Phi = [ b(7)  b(8) ;
            1     0   ];
    
    R   = diag(b(9:14).^2);
    Q = diag([1;zeros(length(Phi)-1,1)]);
    
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
   


