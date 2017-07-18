%=========================================================================
%
%    Simulating the size of the wald test using exponential regression
%    model
%
%=========================================================================
function test_size( )

    clear all;
    clc;
    
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123457) )
          
    beta0  = 1;                                   
    beta1  = 0;                                         
    t      = [5,10,25,100];                                  
    ndraws = 10000;
    wd     = zeros(ndraws,1);
    chi2   = chi2inv(0.95,1);
    output = zeros(2,length(t));
    size   = output(1,:);
    cv     = output(2,:);
      
    for l = 1:length(t)
        j  = 0; 
        x  = randn(t(l),1);     % Explanatory variable (fixed in repeated samples)    

        for i = 1:ndraws           % Main do loop to generate Monte Carlo results   
     
            mue     = exp(beta0 + beta1*x);  %mean
            u       = rand(t(l),1);
            y       = -mue.*log(1 - u);
   
            theta0  = [beta0, beta1];
            options = optimset('LargeScale', 'off', 'Display', 'off');
            theta   = fminunc(@(theta) lnl(theta,y,x),theta0,options);
            H       = Hess(theta,y,x);
   
            %   One restriction
            R       = [0 1];
            Q       = 0;
            wd(i,1) = t(l)*(R*theta' - Q)'*inv( R*inv(-H)*R' )*(R*theta' - Q);

            if wd(i)>chi2
        
                j = j+1;
            end
    
        end
        size(l) = j/ndraws;
        wd_sort = sort(wd,1);
        cv(l) = wd_sort(ndraws*0.95);
    
    end
    disp('                         5        10          25          100 ');
    disp('---------------------------------------------------------------');
    disp(['size:                 ' num2str(size)]);
    disp(['Critical value (5%):  ' num2str(cv)]);

end




%--------------------------- Functions -----------------------------------
% Unconstrained log-likelihood function at each observation
function ln1 = lnlt(theta,y,x)

     beta0 = theta(1);
     beta1 = theta(2);
     mue   = exp(beta0 + beta1*x);
     ln1   = -log(mue) - y./mue;
     
end

% Unconstrained log-likelihood function
function ln = lnl(theta,y,x)

    lnlt1 = lnlt(theta,y,x);
    ln = -mean(lnlt1);

end

% Hessian
function H = Hess(theta,y,x)

beta0  = theta(1);
beta1  = theta(2);
mue    = exp(beta0 + beta1*x);
H      = zeros(2,2);
H(1,1) = -sum(y./mue);
H(1,2) = -sum((x.*y)./mue );
H(2,1) = H(1,2);
H(2,2) = -sum( ((x.^2).*y)./mue );

end