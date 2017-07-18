%=========================================================================
%
%      Program to estimate the consumption capm model by GMM
%
%=========================================================================

function gmm_ccapm( )

    clear all;
    clc;

    % Data are from Hansen and Singleton (adjusted)                   
    load ccapm.mat;
    xdata = ccapm;


    cratio    = xdata(:,1);      % Ratio of real consumption at t+1 and t   
    g_return  = xdata(:,2);      % Real gross return at t+1                        
    r         = g_return - 1;    % Real interest rate at t+1                      
    e         = xdata(:,3);      % Value weighted real returns 
  
    t = length(cratio)-1;           

 
    % Instruments = {const,cratio}   
    %zt = [ ones(t,1)  trimr(cratio,0,1) ]; 
    
    % Instruments = {const,cratio,rlag}
    %zt = [ ones(t,1)  trimr(cratio,0,1) trimr(r,0,1) ];
    
    % Instruments = {const,cratio,rlag,e}
    zt = [ ones(t,1)  trimr(cratio,0,1) trimr(r,0,1) trimr(e,0,1) ];
    
    
    b0 = [0.9;2.0];      
    
    [ bgmm,qmin ] = fminsearch(@(b) q(b,trimr(cratio,1,0),trimr(r,1,0),zt,0),b0);

    % Numerical gradient at optimum
    dg = numgrad(@meaneqn,bgmm,trimr(cratio,1,0),trimr(r,1,0),zt);
       
    % Compute optimal w
    d    = meqn(bgmm,trimr(cratio,1,0),trimr(r,1,0),zt);
    g    = mean(d)';       
    w    = d'*d;
    tau  = 1;
    lmax =0;

    while tau <= lmax
        
        wtau = d((tau+1):size(d,1),:)'*d(1:(size(d,1)-tau),:);
        w    = w + (1.0-tau/(lmax+1))*(wtau + wtau');
        tau  = tau + 1;
        
    end
    w = w/t;
    v = dg'*inv(w)*dg;
    
    % Hansen Sargan J Test
    j_stat = t*q(bgmm,trimr(cratio,1,0),trimr(r,1,0),zt,0);
    

    disp( ['The value of the objective function (q)  = ' num2str(qmin) ] );
    disp( ['J-test is                                = ' num2str(j_stat) ] );
	disp(  'GMM estimates   se     t  ' );
    disp( [ bgmm sqrt(diag(inv(v)/t)) (bgmm./sqrt(diag(inv(v)/t)))] );

end
%
%------------------------- Functions -------------------------------------%
%
%-------------------------------------------------------------------------%
% Define the moment equations 
%-------------------------------------------------------------------------%
function mt = meqn(b,cratio,r,zt)
   
    [~,cols] = size( zt );

    beta = b(1);       %  Relative risk aversion parameter    
    gam  = b(2);       %  Discount parameter                  

    ut = beta*cratio.^(-gam).*(1 + r) - 1;
        
    k = 1;        
    mt = zeros(length(zt),2);
    
    for j = 1:cols      
        mt(:,k) = zt(:,j).*ut;
        k=k+1;
    end             
end
%-------------------------------------------------------------------------%
% Defines the mean of the moment conditions  
%-------------------------------------------------------------------------%
function ret = meaneqn(b,cratio,r,zt)

        ret = mean(meqn(b,cratio,r,zt))';

end
%-------------------------------------------------------------------------%
% GMM objective function which also computes the optimal w   
%-------------------------------------------------------------------------%
function ret = q(b,cratio,r,zt,lmax)

        d = meqn(b,cratio,r,zt);

        g = mean(d)';
        
        w  = d'*d;

        tau=1;

        while tau <= lmax
            wtau = d((tau+1):size(d,1),:)'*d(1:(size(d,1)-tau),:);
            w    = w + (1.0-tau/(lmax+1))*(wtau + wtau');
            tau  = tau + 1;
        end
        t = size(cratio,1);
        w = w./t;
            
        ret = g'*inv(w)*g;
end

