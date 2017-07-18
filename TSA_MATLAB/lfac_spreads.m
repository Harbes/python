%==========================================================================
%
%   Estimate a two factor model of the term structure 
%   using forward spreads relative to the one year yield.
%   Allow for money shocks as an exogenous variable.
%
%==========================================================================
function lfac_spreads( )

    clear all
    clc
    
    % Load the data		
    
    % Yields 1yr, 2yr ... to 30yr    
    data = importdata('daily_finance.csv');
    
    yields  = data(:,1:30);
    forwards = data(:,31:end-1);
    money    = data(:,end);
     
    % Compute spreads relative to the 1yr yield of the (n-1) forwards maturing at n				
    y = log(1+forwards(:,[1 3 5 7 9])/100) - repmat(log(1+yields(:,1)/100),1,5);      
    y = 10000*y;                                                       
    x = money;                                                           
    y = bsxfun(@minus,y,mean(y));
    t = length(y);

    % Estimate the model 
    flag = 1;
    start = [
                1.3767001935080581 ;
                6.2447235611066363 ;
                8.3966862784726608 ;
                9.2782017956204559 ;
                9.5673550480256875 ;
                2.1681182372527874 ;
                2.4537907722377326 ;
                0.4868549700052381 ;
                -1.2716822171291355 ;
                -2.5123991860804984 ;
                11.2849486658020870 ;
                -0.7663744607435562 ;
                 5.0545819842149111 ;
                 1.1823456636092726 ;
                -8.3990683439717611 ;
                 4.1675429449096004 ;
                 2.7998959843657421 ;
                -3.8632689276557515 ;
                 0.6123722675895771 ];
    
     ops   = optimset('LargeScale','off','Display','iter', ...
                    'MaxFunEvals',Inf,'MaxIter',Inf);
 
    [bhat] = fminunc(@(b) neglog(b,y,x,flag),start,ops);
    
    % Restimate without restrictions to get standard errors
    flag  = 0;
    start = [bhat(1:15); tanh(bhat(16:17)); bhat(18:19)];
        
    [ bhat,lf,~,~,~,hess ] = fminunc(@(b) neglog( b,y,x,flag ),start,ops );
    
    vc = (1/t)*inv(hess);
    
    disp( ['Log-likelihood function = ',num2str(lf) ] );
    disp( ' ' )   
    disp( 'Params    Std. Errors ' )
    disp( [bhat  sqrt(diag(vc))] )

end

%--------------------------- Functions ----------------------------------
% 
%--------------------------------------------------------------------------
%   Wrapper function to set up and call the Kalman filter
%--------------------------------------------------------------------------
function lf = neglog(b,y,x,flag)

    % Unpack the parameter vector
    Lam = [ b(1:5) , b(6:10)];
    Gam = b(18:19);
    if flag
        Phi = diag(tanh(b(16:17)));
    else
        Phi = diag(b(16:17));
    end

    %R = eye(6);
    R = diag(b(11:15).^2);
    Q = eye(2);
        
    lf = -mean( kalman(y,x,Phi,Lam,Gam,R,Q) );
   
end
%--------------------------------------------------------------------------
% Kalman filter
%--------------------------------------------------------------------------
function lnl = kalman(y,x,Phi,Lam,Gam,R,Q)


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
        st = Phi*s0 + Gam*x(i);                    
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



