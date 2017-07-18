%=========================================================================
%
%		Finite sample properties of the Bai and Ng (2004) of the PANIC model
%
%=========================================================================
function lfac_panic( )

    clear all
    clc
    
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) );

    
    % Settings
    t = 200;                % Number of observations                   
    n = 50;                 % Number of dependent variables        
    k = 1;                  % Number of factors     
    ndraws = 5000;          % Number of draws   
    lam = 1 + randn(n,1);   % Choose lambda parameters N(1,1)          
    delta = 1.0;            % Choose delta parameter (same for all u) 

    % Allocate arrays
    lam_sw = zeros(ndraws,n);  % Loadings of the Stock-Watson estimator
    lam_bn = zeros(ndraws,n);  % Loadings of the Bai-Ng estimator


    for iter = 1:ndraws
        % Simulate the data       
        w  = randn(t,n);
        u  = recserar(w,zeros(1,n),delta*ones(1,n));  % Idiosyncratics                 
        v  = randn(t,1);
        s  = recserar(v,0,1);                         % Factor                              
        s  = (s - mean(s))./std(s);                   % Standardize factor                   
        s0 = s;                                       % Define true factor 
        y  = bsxfun(@times,s,lam') + u;
         
        % First difference the data     
        dy = trimr(y,1,0) - trimr(y,0,1);
        dy = bsxfun(@minus,dy,mean(dy));
                
        % Principal component analysis using levels  (Stock Watson) 
        [ ~,s ] = princomp( y );                
        s        = s(:,1:k);                % First k principal components
        s        = (s - mean(s))./std(s);
                
        % Principal component analysis using first differences  (Bai-Ng estimator)
        [ ~,ds ] = princomp( dy );                
        ds       = ds(:,1:k);                % First k principal components
        ds       = (ds - mean(ds))./std(ds);
 
        % Estimate factor loadings  
        b = ([ones(length(y),1) s])\y;  
        %  Stock-Waton: using levels (pick out the slopes)  
        lam_sw(iter,:) = b(2,:);     

        % Bai-Ng: using first differences
        tmp            = (ds\dy)';
        lam_bn(iter,:) = tmp(:);         
        
        
    end
    % Some plots
    figure(1)
    plot(1:length(s),[s0 s]);
    title('Stock-Watson')
 
    figure(2)
    plot(1:length(s), [s0 [s0(1); cumsum(ds)] ]);
	title('Bai-Ng')
  
    % Compute statistics of sampling distributions and print results     
    mse_sw = mean( bsxfun(@minus,lam_sw,lam').^2 );
    mse_bn = mean( bsxfun(@minus,lam_bn,lam').^2 );

    disp(' ')
    disp(['T       = ', num2str(t) ]);
    disp(['delta   = ', num2str(delta) ]);
    disp(' ');

    disp('Stock-Watson Estimator');
    disp('Mean   Bias   MSE '),  
    disp( [mean(lam_sw)' mean(lam_sw)'-lam  mse_sw' ] );
   

    disp('Bai-Ng Estimator');
    disp('Mean   Bias   MSE '),  
    disp( [mean(lam_bn)' mean(lam_bn)'-lam  mse_bn' ] );

    disp(' ')
    disp('Stock-Watson Estimator Overall');
    disp(['Mean    = ' num2str(mean(mean(lam_sw))) ]);
    disp(['Bias    = ' num2str( mean(mean(lam_sw)'-lam) ) ]); 
    disp(['MSE     = ' num2str( mean(mse_sw) ) ]);

    disp(' ')
    disp('Bai-Ng Estimator Overall');
    disp(['Mean    = ' num2str(mean(mean(lam_bn))) ]);
    disp(['Bias    = ' num2str( mean(mean(lam_bn)'-lam) ) ]); 
    disp(['MSE     = ' num2str( mean(mse_bn) ) ]);

end
