% ========================================================================
%
%      Monte Carlo properties of GMM estimator of the gamma distribution
%
% ========================================================================
function gmm_gammasim( )

    clear all
    clc

    format short
    
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) );

    % Set the parameters of the Monte Carlo experiment  
    reps = 200;
    
    % Now uncomment only one of the following options
    % Estimation  
    %dist = 'gam'; t = 50; a0v = [1,2,3,4,5];
    dist = 'gam'; t = 100; a0v = [1,2,3,4,5];
    % //dist = 'gam'; t = 200; a0v = [1,2,3,4,5];
    % //dist = 'gam'; t = 400; a0v = [1,2,3,4,5];
    
    % Testing: size
    %dist = 'gam'; t = 100; a0v = 1;
    %dist = 'gam'; t = 100; a0v = 2;
    %dist = 'gam'; t = 100; a0v = 3;
    %dist = 'gam'; t = 100; a0v = 4;
    %dist = 'gam'; t = 100; a0v = 5;
    
    %dist = 'gam'; t = 200; a0v = 1;
    %dist = 'gam'; t = 200; a0v = 2;
    %dist = 'gam'; t = 200; a0v = 3;
    %dist = 'gam'; t = 200; a0v = 4;
    %dist = 'gam'; t = 200; a0v = 5;
    % 
    %Testing: power
    %dist = 'gam'; t = 200; a0v = [1,1.05,1.1,1.15,1.2,1.25,1.3];
     
    % Misspecification
    %dist = 'exp'; t = 200; a0v = [1,1.2,1.4,1.6,1.8,2.0];
    %dist = 'exp'; t = 400; a0v = [1,1.2,1.4,1.6,1.8,2.0];
    %dist = 'exp'; t = 800; a0v = [1,1.2,1.4,1.6,1.8,2.0];
     
    %dist = 'exp'; t = 200; a0v = 1;
    %dist = 'exp'; t = 200; a0v = 1.2;
    %dist = "exp"; t = 200; a0v = 1.4;
    %dist = 'exp'; t = 200; a0v = 1.6;
    %dist = 'exp'; t = 200; a0v = 1.8;
    %dist = 'exp'; t = 200; a0v = 2;
   
    simexp(t,reps,dist,a0v);
    
end
    
%
%--------------------------- Functions -----------------------------------
%    
%-------------------------------------------------------------------------
%  Monte Carlo experiment for different parameter inputs
%-------------------------------------------------------------------------

function simexp(t,reps,dist,a0v)

    cols = length(a0v);
    
    % Allocate arrays
    aML  = zeros(reps,cols); 
    a1   = zeros(reps,cols);  
    a2   = zeros(reps,cols);
    bML  = zeros(reps,cols);
    b1   = zeros(reps,cols);
    b2   = zeros(reps,cols);
    seML = zeros(reps,cols); 
    se1  = zeros(reps,cols);  
    se2  = zeros(reps,cols);
    seQ2 = zeros(reps,cols);
    tML  = zeros(reps,cols); 
    t1   = zeros(reps,cols);  
    t2   = zeros(reps,cols);
    tQ2  = zeros(reps,cols);
    j2   = zeros(reps,cols);
    
    for ac = 1:length(a0v)

        for rep = 1:reps    
            if strcmp(dist,'gam')
                y = gamrnd(a0v(ac),1,t,1);
            elseif strcmp(dist,'exp')
                y = -log(rand(t,1))*a0v(ac);
            else
                disp('Wrong distribution')
            end
        
            % Maximum likelihood
            ops = optimset( 'LargeScale', 'off','Display', 'off' );
            [aML(rep,ac),~,~,~,~,H]=fminunc( @(a) neglog(a,y),mean(y),ops );
    
            bML(rep,ac)  = aML(rep,ac)-a0v(ac);
            seML(rep,ac) = sqrt(1/(t*H));
            tML(rep,ac)  = bML(rep,ac)/seML(rep,ac);
            
            % GMM using the first moment
            a1(rep,ac)  = mean(y);           
            b1(rep,ac)  = a1(rep,ac)-a0v(ac);
            se1(rep,ac) = std(y)/sqrt(t);
            t1(rep,ac)  = b1(rep,ac)/se1(rep,ac);
            
            % GMM using two moments
            [a2(rep,ac),qmin(rep,ac),~,~,~,H]=fminunc( @(a) gmmcrit(a,y),mean(y),ops );

            b2(rep,ac)  = a2(rep,ac)-a0v(ac);
            se2(rep,ac) = sqrt(1/(t*H));
            t2(rep,ac)  = b2(rep,ac)/se2(rep,ac);
            
            m            = [(y-a2(rep,ac))  (y.^2-a2(rep,ac)*(a2(rep,ac)+1)) ];
            WTi          = inv(m'*m/length(m));
            D            = [ 1; (2*a2(rep,ac)+1) ];
            seQ2(rep,ac) = 1/sqrt(t*D'*WTi*D);
            tQ2(rep,ac)  = b2(rep,ac)/seQ2(rep,ac);
            
            j2(rep,ac)   = 2*t*qmin(rep,ac);
 
        end
        
    end
    disp([' Distribution    ', dist] ) 
    disp('Results for Maximum Likelihood')
    disp('    Bias     Std. err.  t >1.96 ')
    disp([mean(bML)' mean(seML)' mean(abs(tML)>1.96)' ])
    disp(' ')
    disp('Results for GMM (1 moment condition)')
    disp('Results for Maximum Likelihood')
    disp('    Bias     Std. err.  t >1.96 ')
    disp([mean(b1)' mean(se1)' mean(abs(t1)>1.96)' ])
    disp(' ')
    disp('Results for GMM (2 moment conditions)')
    disp('Results for Maximum Likelihood')
    disp('    Bias     Std. err.  t>1.96  Std.err(Q)  t(Q)>1.96 J>3.84')
    disp([mean(b2)' mean(se2)' mean(abs(t2)>1.96)' mean(seQ2)' mean(abs(tQ2)>1.96)' mean(j2>3.84)'])

end
    
%-------------------------------------------------------------------------
% Negative log-likelihood function for gamma distribution with beta = 1
%-------------------------------------------------------------------------

function logl = neglog( a,y )


    logl =  (a-1)*mean(log(y)) - mean(y) - gammaln(a);
    logl = - logl;

end
%-------------------------------------------------------------------------
% GMM criterion function for gamma distribution with beta = 1
%-------------------------------------------------------------------------
function q = gmmcrit(a,y)

      m = [ (y-a)  (y.^2 - a*(a+1)) ];
      w = m'*m/length(y);
      q = 0.5*mean(m)*inv(w)*mean(m)';

end

