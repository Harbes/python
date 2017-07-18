%=========================================================================
%
%   Approximate a nonlinear time series model using a 
%   variable parameter model (Granger (2008))
%
%=========================================================================
function nlm_tvarying( )

    clear all;
    clc;

    % Initialise the random number generator  
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',12) ); 	

    % Simulate data	
    t = 200;
    y = zeros(t+100,1) ;
    v = sqrt(3)*randn(t+100,1);

    % Parameters
    % Choose gam = 0.0 for first experiment and gam = 0.8 for the second experiment
    phi0 = 10.0;
    phi1 = 0.6;
    phi2 = 0.3;
    gam  = 0.8;            


	for k = 3:t+100

    y(k) = phi0 + phi1*y(k-1) + phi2*y(k-2) + gam*y(k-1)*v(t-1) + v(k);
    
    end
    y = trimr(y,100,0);
    
    % True conditional expectation
    my_true = phi0 + phi1*trimr(y,1,1) + phi2*trimr(y,0,2);                 

    % Nonparametric approximation     
    yt = trimr(y,1,0);
    x  = trimr(y,0,1);
    xt = x;
    h  = 1.06*std(yt)'*t^(-1/5);

    fx  = mean( normpdf( (repmat(x,1,length(x)) - repmat(xt',length(x),1))/h )'/h )';
   	fyx = mean( normpdf( (repmat(x,1,length(x)) - repmat(xt',length(x),1))/h )'.*repmat(yt,1,length(yt)) )'/h;

    my_npr = fyx./fx;
    my_npr = trimr(my_npr,1,0);

    % Estimate the model as a Kalman filter with time-varying parameters      
    ops     = optimset('LargeScale','off','Display','iter');
    theta_0 = [0.1 ; 0.1];
    smooth  = 0;             % Toggle to compute the smoothed factor (smooth = 1).
    theta = fminunc(@(b) neglog(b,y,smooth),theta_0,ops);


    % Compute the smoothed factor estimates using the constrained MLEs  
    smooth = 1;
    [~,s_smooth] = loglt(theta,y,smooth);

    my_tvar = trimr(s_smooth,1,0).*trimr(y,0,1);                   
    my_tvar = trimr(my_tvar,1,0);

    % Graph conditional expectations
    figure(1);
    clf;
    tt = 1:length(my_true);
    plot(tt,my_true,'r',tt,my_tvar,'g',tt,my_npr,'b');
    ylabel('Conditional mean');
    xlabel('y');
    legend('True','Time Varying','Nonparametric')

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Negative log-likelihood function
%-------------------------------------------------------------------------
function lf = neglog(b,y,smooth)

    [tmp,~] = loglt(b,y,smooth);
    lf = -mean( tmp );
    
end
%-------------------------------------------------------------------------
%  Kalman filter
%-------------------------------------------------------------------------

function [logl,s_smooth] =  loglt(b,y,smooth)

    [t,n]= size(y);
    logl = zeros(t,1);

    % Define the matrices for the filter   
    lam = 0.0;
    r   = diag(b(1).^2);
    phi = tanh(b(2));
    q   = 1; 
    k  = size(q,1);

      s_predict = zeros(t,k);     %     Declare predicted factor estimates     
      s_update  = zeros(t,k);     %     Declare updated factor estimates       
      s_smooth  = zeros(t,k);     %     Declare smoothed factor estimates      
      p_predict = zeros(t,k^2);   %     Declare updated factor estimates vcov  
      p_update  = zeros(t,k^2);   %     Declare smoothed factor estimates vcov 
      p_smooth  = zeros(t,k^2);   %     Declare smoothed factor estimates vcov 


    % Run through the filter and construct the likelihood 
    st = zeros(k,1);
    pt = eye(k)*0.1;    %     Initialization based on the diffuse prior distribution   

    s0 = st;
    p0 = pt;

    for i = 2:t
    
         % Prediction      
         st = phi*s0;
         pt = phi*p0*phi' + q;

         % Observation     
         lam = y(i-1,:);     
         mt  = lam*st;
         vt  = lam*pt*lam' + r;
         ut  = y(i,:)' - mt;

         % Updating        
         s0 = st + pt*lam'*inv(vt)*ut;
         p0 = pt - pt*lam'*inv(vt)*lam*pt;

         % Log-likelihood  
         logl(i) = - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5*ut'*inv(vt)*ut;

         if smooth == 1

            s_predict(i,:) = st';       %     Store predicted factor estimates        
	        s_update(i,:)  = s0';       %     Store updated factor estimates          
	        p_predict(i,:) = vec(pt)';  %     Store predicted factor estimates vcov   
	        p_update(i,:)  = vec(p0)';  %     Store updated factor estimates vcov     

         end
    end
    
    % Generate smoothed factor estimates
    if smooth == 1
        s_smooth = s_update;
        p_smooth = p_update;

        for i = 1:t-1
            j = reshape(p_update(t-i,:),k,k)'*phi'*inv(reshape(p_predict(t-i+1,:),k,k)'); 
            s_smooth(t-i,:) = ( s_update(t-i,:)' + j*(s_smooth(t-i+1,:)' - s_predict(t-i+1,:)') )';
            p_smooth(t-i,:) = vec( reshape(p_update(t-i,:),k,k)' + j*(reshape(p_smooth(t-i+1,:),k,k)' - reshape(p_predict(t-i+1,:),k,k)')*j' )';
        end      
    end
end

      




         


