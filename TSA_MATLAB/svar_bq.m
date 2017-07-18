%==========================================================================
%
%      Estimate the Blanchard-Quah SVAR model
%
%==========================================================================
function svar_bq( )

    clear all
    clc
    
    % Load data (1950Q1 - 2009Q3)
    %    (1) Real GDP, s.a.
    %    (2) Unemployment rate, s.a.
    load('bq_us.mat');
    %load('bq_uk.mat');
    %load('bq_jp.mat');

    ytdata = ytdata_us;
    
    % Construct variables for use in VAR
    lrgdp  = log(ytdata(:,1));
    urate = ytdata(:,2);

    y = [400*(trimr(lrgdp,1,0) - trimr(lrgdp,0,1))  trimr(urate,1,0) ];
    [t,n] = size(y);
    
    disp('Covariance matrix of the data')
    disp( cov(y) );

    % Estimate the VAR with p lags, a constant and a time trend		**/
    p = 8;		% Order of VAR lags	 
    q = 40;		% Order of VMA lags		 

    ylag = lagmatrix(y,1:p);
    ylag = [ones(t,1) ylag];
    ylag(any(isnan(ylag),2),:) = [];
    
    yt  = trimr(y,p,0);
    bar = ylag\yt;
    v   = yt - ylag*bar;
    vc  = v'*v/length(v);
    
    disp('Covariance matrix of VAR residuals')
    disp( vc );

    % Constuct A(L) (VAR) matrices and hence the C(1) long-run matrix
    bar = trimr(bar,1,0);
    k   = size(v,2);
   	a   = zeros(k^2,p);
    a1  = eye(k);

    for i = 1:p
        
        tmp    = bar(1+k*(i-1):i*k,:);
    	a(:,i) = tmp(:);
	    a1     = a1 - reshapeg(a(:,i),k,k);
    end
      
    % Invert A(1) matrix needed for the MLE
    a1inv = inv(a1); 
    disp('Inverse of a1 matrix')    
    disp(a1inv);

    % Esimate Blanchard-Quah SVAR model with long-run restrictions
    ops = optimset(    'LargeScale',           'off', ...
                        'Display',              'off', ...
                        'MaxIter',              2000,  ...
                        'MaxFunEvals',          4000 );        
    bstart     = vech(chol(vc)');         
    [ b,fval,aa,aa,aa,H]  = fminunc( @(b) neglog(b,v,a1),bstart,ops );     
 
    disp(['Log-likelihood function  = ',num2str(-fval) ]);
    
    vcov = (1/length(v))*inv(H);
    
    % Test of Okun's Law: long-run version 
    c = b(2)/b(1) + 3;              % Aggregate supply shock parameters 
    d = [-b(2)/b(1)^2   1/b(1)  0];
    w = c'*inv(d*vcov*d')*c;

    disp('Wald test of Okuns Law (aggregate supply version)');
    disp(['Mean            = ',num2str(b(2)/b(1)) ]);
    disp(['Wald statistic  = ',num2str(w) ]);
    disp(['p-value         = ',num2str(1-chi2cdf(w,1)) ]);
    disp(' ' );
    disp('Wald test of Okuns Law (aggregate demand version)');
    disp('Aggregate demand shocks have no affect on output in the long-run');
    disp(' ');
    
    % Impuse responses
    % Construct C(L) matrices (vector moving average) 
    c = eye(k);
    c = c(:);

    for i = 1:q

       ss = zeros(k,k);
       j   = 1.0;
        
       while j <= min( [ p i ]);

          ss = ss + reshapeg(a(:,j),k,k)*reshapeg(c(:,i-j+1),k,k);
          j   = j + 1;
       end
       tmp = ss';
       c  = [ c  tmp(:) ];

    end
    
    
    % Construct orthogonal impulse response functions
    s = svarmat(b,a1);
    impulse = vec(s');

    for i=2:q+1

       tmp      = ( reshapeg( c(:,i),k,k )*s )';
       impulse  = [ impulse tmp(:) ];

    end 

	impulse = impulse';
    
    % Compute cumulative sum of impulse responses for levels 
	impulse = [ cumsum(impulse(:,[1 2])) impulse(:,[3 4]) ];		

    % Construct variance decompositions 
    tmp0   = reshapeg( cumsum( impulse(1:end-1,:).^2 ),q*k,k );
    tmp1   = repmat( sum(tmp0,2),1,k );
    decomp = reshapeg( 100*tmp0 ./ tmp1 , q , k^2 );
    
    % Plot impulse responses 
    figure(1);
    clf
   
    subplot(2,2,1)
    plot(seqa(0,1,q+1),[impulse(1:q+1,1) zeros(q+1,1)] )
	title('Agg. supply shock');
	ylabel('Real Output');
    xlabel ('Quarter');

    subplot(2,2,2)
    plot(seqa(0,1,q+1),[impulse(1:q+1,2) zeros(q+1,1)] );
	title('Agg. demand shock');
	ylabel('Real Output');
    xlabel ('Quarter');

    subplot(2,2,3)
    plot(seqa(0,1,q+1),[impulse(1:q+1,3) zeros(q+1,1)] );     
	title('Agg. supply shock');
	ylabel('Unemployment');
    xlabel ('Quarter');
 
    subplot(2,2,4)
	plot(seqa(0,1,q+1),[impulse(1:q+1,4) zeros(q+1,1)]);   
	title('Agg. demand shock');
	ylabel('Unemployment');
    xlabel ('Quarter');

          
    % SVAR with long-run restriction and Okuns Law restriction   
	bstart = b([1 3]);
    [ b,fval,aa,aa,aa,H]  = fminunc( @(b) neglogc(b,v,a1),bstart,ops );     

    disp(['Log-likelihood function  = ',num2str(-fval) ]);
    disp(' ');
    

    
end
%--------------------------- Functions ----------------------------------
% 
%--------------------------------------------------------------------------
% Log-likelihood function for SVAR with long-run restrictions
%--------------------------------------------------------------------------
function lf = neglog(b,v,a1)

    [t,n] = size(v);
    lnl   = zeros(t,1);

    f   =   [ b(1)        0   ;
              b(2)     -abs(b(3))  ];      

	s   = a1*f;		% Long-run restrictions		 
	vc  = s*s';
    ivc = inv(vc);
    
	for i=1:t
		lnl(i) = -0.5*n*log(2*pi) - 0.5*log(det(vc)) ...
                 - 0.5*v(i,:)*ivc*v(i,:)';    
    end
    lf = -mean( lnl );

end
%--------------------------------------------------------------------------
% Return SVAR matrices
%--------------------------------------------------------------------------
function s = svarmat(b,a1)
    
    f   =   [ b(1)        0   ;
              b(2)     -abs(b(3))  ];      

	s   = a1*f;		% Long-run restrictions		 
end


%--------------------------------------------------------------------------
% Log-likelihood function for SVAR with Okuns Law and long-run restrictions
%--------------------------------------------------------------------------
function lf = neglogc(b,v,a1)

    [t,n] = size(v);
    lnl   = zeros(t,1);

    f   =   [  b(1)        0   ;
             -3*b(1)     -abs(b(2))  ];      

	s   = a1*f;		% Long-run restrictions		 
	vc  = s*s';
    ivc = inv(vc);
    
	for i=1:t
		lnl(i) = -0.5*n*log(2*pi) - 0.5*log(det(vc)) ...
                 - 0.5*v(i,:)*ivc*v(i,:)';    
    end
    lf = -mean( lnl );

end

