%=========================================================================
%
%   Estimate a MIDAS regression model for the hedge funds 
%
%=========================================================================
function garch_midas( )

    clear all
    clc
    
    % Load daily data: March 31st 2003 to Aug. 4th 2011 
    load daily_hedge;
    
    % Compute daily returns
    rhd = trimr(log(p_hedge),1,0) - trimr(log(p_hedge),0,1);                  
    rmd = trimr(log(p_market),1,0) - trimr(log(p_market),0,1);                 

    % Compute monthly returns from daily returns 
    % assuming a month has 22 days (starting with April 2003)
    rh = sum( reshapeg( rhd,length(rhd)/22,22),2 );
    rm = sum( reshapeg( rmd,length(rmd)/22,22),2 );

    %  Lag length to compute the MIDAS 
    % lag=1 -> conditional variance based on daily returns in the previous month
    % lag=2 -> conidtional variance based on daily returns in the previous two months         
    lag = 18;           
    tm = length(rh);           %  Number of months in the monthly data set
    td = length(rhd);          %  Number of days in the daily data set
    n = round(td/tm);          %  Number of days in a month

    % Squared daily returns scaled by 22 to convert to monthly returns
    % already lagged given the way that the daily data are defined 
    rhd2 = 22*rhd.^2;                  


    % Construct global financial crisis dummy variable
    d_gfc          = zeros(tm,1);    
    [ ~,k ]  = min(rh);
    d_gfc(k) = 1;

    % Estimate the MIDAS model
    ops   = optimset('LargeScale','off','Display','iter');
    start = 0.1*ones(6,1);
    
    [theta,lf1,~,~,~,hess] = fminunc(@(b) neglog(b,rh,rm,rhd2,d_gfc,tm,n,lag),start,ops);
   
    clear start
    lf1 = -lf1;
    vc  = (1/(tm-lag))*inv(hess);
    
    % Plot estimated weight function 
    w = getweights(theta(5),theta(6),n*lag);   
    plot(w)
    ylabel('MIDAS weights');
    xlabel('i');

    
    % Wald test of theta(5)=theta(6)=0
    r = [ 0  0  0  0  1  0 ;
          0  0  0  0  0  1 ];
    q = [ 0 ; 0 ];
    w = (r*theta - q)'*inv(r*vc*r')*(r*theta - q);

    disp(['Wald statistic     = ',num2str(w)]);
    disp(['p-value            = ',num2str(1-chi2cdf(w,2))]);
 
end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Log-likelihood function
%-------------------------------------------------------------------------
function lf = neglog(b,rh,rm,rhd2,d_gfc,tm,n,lag)
    
    theta1 = b(5);
    theta2 = b(6);    
        
    % Precompute weights
    maxlag = n*lag;
    w      = getweights(theta1,theta2,maxlag);
    
    % Now need 3 counters, k,i,j to compute midas variance
    % k is a simple placeholder for the calculated monthly variance
    % i loops over the daily data and is incremented by the number of days in month
    % j loops from 1 to maxlag over the daily data to compute the required sum  
    
    vm = zeros(tm,1);
    
    % If maxlag is greater than number of days in a month then pad the returns
    if maxlag > n 
        rhd2 = [ zeros(maxlag-n,1); rhd2];
    end    
    
    k = 1;
    for j = maxlag+1:n:length(rhd2)+1      
        for i = 1:maxlag       
            vm(k) = vm(k) + rhd2(j-i)*w(i); 
        end
        k=k+1;
    end       
    % Trim extraneous start up months and lag the result  
    if lag > 1
        vtmp = trimr(vm,lag-1,1);
    else
        vtmp = trimr(vm,0,1);
    end
        
    v  = trimr(rh,lag,0) - b(1) - b(2)*trimr(rm,lag,0) ...
        - b(3)*vtmp - b(4)*trimr(d_gfc,lag,0);   
    s2 = v'*v/length(v);                                              
    z  = v./sqrt(s2);
    f  = - 0.5*log(2*pi) - 0.5*log(s2) - 0.5*z.^2;   
    lf = -mean( f );

end

%-------------------------------------------------------------------------
%  Return weights normalised to sum to unity
%-------------------------------------------------------------------------
function w = getweights(theta1,theta2,maxlag)

    i = 0:1:maxlag-1;
    %i = seqa(0,1,maxlag);
    w = exp( theta1*i/1.0e2 + theta2*i.^2/1.0e4  );
    w = w./sum(w);    
    
end
