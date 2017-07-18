%=========================================================================
%
%      Lee, Granger, White neural network test J Ects (1993) 
%
%=========================================================================
function nlm_lgwtest( )

    clear all;
    clc;

    % Initialise the random number generator  
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',12345) ); 	

    % Parameters
    reps = 2000;       % Number of draws                 
    q    = 2;          % Number of activation functions    
    nobs = 250; 


    %  Null hypothesis: AR Model   
    ntest = zeros(reps,1);
  
    for r = 1:reps

        y = zeros(nobs+100,1);
        v = randn(nobs+100,1);
        for t=2:nobs+100;
            y(t) = 0.6*y(t-1) + v(t);          
        end
        y = trimr(y,100,0);   
        x = trimr(y,0,1);    
        y = trimr(y,1,0);   


        [ypred,ntest(r),pv] = neural(y,[ones(size(x,1),1) x],q);

    end
    
    dist = quantile(ntest,[.010 .025 .050 .100 .500 .900 .950 .975 .990]);

    disp(' ')
	disp('Empirical distribution (under the null of linearity)');
 
	disp(['1.0 per cent      =  ', num2str(dist(1)) ]);
    disp(['2.5 per cent      =  ', num2str(dist(2)) ]);
    disp(['5.0 per cent      =  ', num2str(dist(3)) ]);
    disp(['10.0 per cent     =  ', num2str(dist(4)) ]);
    disp(['50.0 per cent     =  ', num2str(dist(5)) ]);
    disp(['90.0 per cent     =  ', num2str(dist(6)) ]);
    disp(['95.0 per cent     =  ', num2str(dist(7)) ]);
    disp(['97.5 per cent     =  ', num2str(dist(8)) ]);
    disp(['99.0 per cent     =  ', num2str(dist(9)) ]);

    cv = dist(7);     

	disp(' ' );


    %  Alternative hypothesis: Bilinear Model   
    ntest = zeros(reps,1);
  
    for r = 1:reps

        y = zeros(nobs+100,1);
        v = randn(nobs+100,1);

        for t=3:nobs+100;

            y(t) = 0.7*y(t-1)*v(t-2) + v(t);          
        end

        y = trimr(y,100,0);   


        x = trimr(y,0,1);    
        y = trimr(y,1,0);   


        [ypred,ntest(r),pv] = neural(y,[ones(size(x,1),1) x],q);

    end
    
    dist = quantile(ntest,[.010 .025 .050 .100 .500 .900 .950 .975 .990]);

    disp(' ')
	disp('Empirical distribution (under the alternative)');
 
	disp(['1.0 per cent      =  ', num2str(dist(1)) ]);
    disp(['2.5 per cent      =  ', num2str(dist(2)) ]);
    disp(['5.0 per cent      =  ', num2str(dist(3)) ]);
    disp(['10.0 per cent     =  ', num2str(dist(4)) ]);
    disp(['50.0 per cent     =  ', num2str(dist(5)) ]);
    disp(['90.0 per cent     =  ', num2str(dist(6)) ]);
    disp(['95.0 per cent     =  ', num2str(dist(7)) ]);
    disp(['97.5 per cent     =  ', num2str(dist(8)) ]);
    disp(['99.0 per cent     =  ', num2str(dist(9)) ]);

    cv = dist(7);     

	disp(' ' );

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Estimate and test a neural network
%-------------------------------------------------------------------------
function  [ pred,ntest,pval ] = neural(y,x,q)

    activate = zeros(size(y,1),q);

    for i = 1:q;

       gam    = 4*rand(size(x,2),1)-2.0;  %  Uniform random numbers [-2,2]  )
       lambda = x*gam;
       
       activate(:,i) = (1+exp(-lambda)).^(-1);
       
    end
    xx    = [x activate];
    uhat  = y - x*(x\y);                        
    ehat  = uhat - xx*(xx\uhat);  
    ntest = size(y,1)*(1 - sum(ehat.^2)/sum(uhat.^2));
    pred  = xx*(xx\y);         
    pval  = 1 - chi2cdf(ntest,q);

end

