%==========================================================================
%
%   	Program to estimate a Structural VAR for the share example 
%       using Blanchard-Quah long-run restrictions
%		with present value long-run restriction.
%
%==========================================================================

function svar_port( )

    clear all
    clc
    
    % Load data
        % Data from 1980:1 to 2002:2, taken from the RBA Bulletin database    
		%1. US Share price index, average over the quarter, S&P500
		%2. US/AU exchange rate, average over the quarter
		%3. CPI, all groups, 1989-1990 =100, seasonally adjusted
		%4. Nominal GDP, expenditure, total, s.a.
		%5. Real GDP, expenditure, total, 2000/20001 $m, s.a.
		%6. Australian share price index, average over the quarter, S&P/ASX 200
        %7. 90 day bank accepted bill interest rate, average over the quarter


    load portfolio_data;

    % Data transformations
    ipd   = ytdata(:,4)./ytdata(:,5);		% Implicit price deflator
    r3    = ytdata(:,7)/100;
    lipd  = log(ipd);
    lcpi  = log(ytdata(:,3));
    lrao  = log(ytdata(:,6)./ipd);
    lr3   = log(r3);
    lrgdp = log(ytdata(:,5));
    lrdj  = log(ytdata(:,1)./(ytdata(:,2).*ipd));	% US share price converted into Australian dollars and then expressed in reals 
    inf   = 1*(trimr(lipd,4,0) - trimr(lipd,0,4));

    yact  = 100*( [ lrgdp lr3 lrao lcpi lrdj ] );
    y     = trimr(yact,1,0) - trimr(yact,0,1);

    % Create dummy variables
    d87        = zeros(102,1);	
    d87(32)    = 1.0;			% Stock market crash dummy variable
    d2000      = zeros(102,1); 
    d2000(83)  = 1.0;           % GST introduction in 2000:3c

    % SVAR Parameters
    p = 2;		% VAR order
    q = 40;		% VMA order


    % Estimate the VAR with p lags
  	ylag     = [ ones(length(y),1) trimr( [d87 d2000],1,0) ];
    nconst   = size( ylag,2 );

    for i = 1:p

        ylag = [ trimr(ylag,1,0) trimr(y,0,i) ];

    end

    bar   = ylag\trimr(y,p,0);
   	mu    = bar(1:nconst,:);
   	v     = trimr(y,p,0) - ylag*bar;
    omega = v'*v/size(v,1);
    
    % Constuct A(L) (VAR) matrices and hence the C(1) long-run matrix
    bar = trimr(bar,nconst,0);
    k   = size(v,2);
   	a   = zeros(k^2,p);
    a1  = eye(k);

    for i = 1:p
        
        tmp    = bar(1+k*(i-1):i*k,:);
    	a(:,i) = tmp(:);
	    a1     = a1 - reshapeg(a(:,i),k,k);
    end

    % Estimate unrestricted SVAR
      options = optimset(   'LargeScale',           'off', ...
                            'Display',              'off', ...
                            'MaxIter',              1000,   ...
                            'MaxFunEvals',          4000 ); 
  
    bstart = [
                1.0389699473962597 
                4.3007505638494345 
                5.1966616592138655 
                10.0996788841097480 
                2.7076795082063354 
                1.4379097951392796 
                -0.7123778562836982 
                0.7750813995549357 
                1.6190786936559989 
                1.3494397321299636 
                10.7524502796372890 
                3.5142132747486432 
                -5.1782573061530375  ];
    [ b,fvalu]  = fminunc( @(b) neglogu(b,v,a1),bstart,options );     

    fvalu = -fvalu
    disp(' Unrestricted Log-likelihood function ');
    disp( fvalu );


    % Estimate restritced SVAR by maximum likelihood 
    bstart = [
                1.038918123051897 
                4.301175420613745 
                5.180760075404540 
                10.10913871549879 
                2.707544494664958 
                1.428884030946225 
               -0.7120443107363003 
                0.7725526070345229 
                1.620606390960305 
                1.349442082771462 
                10.75208272049359 
                3.513290773831895            ];
                                              
    [ b,fvalr]  = fminunc( @(b) neglogl(b,v,a1),bstart,options );     
   
    fvalr = -fvalr;
    disp(' Restricted Log-likelihood function ');
    disp( fvalr );
    
    % Likelihood ratio test
    lr = -2*length(v)*(fvalr-fvalu);
    
    disp(['LR test        = ',num2str(lr)] );
    disp(['p-value        = ',num2str(1-chi2cdf(lr,1)) ]);
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
    % Long run restrictions at optimal parameters
    f   =     [   b(1)   0       0       0      0       ;
        	      b(2)   b(3)    b(4)    0      0       ;
	              b(5)  -b(3)    b(6)    0      b(12)   ;
   		          b(7)   b(8)    b(9)    b(10)  0       ;
			       0     0       0       0      b(11) ] ;

	s      = a1*f;		% Long-run restrictions
    
    impulse = vec(s');

    for i=2:q+1

       tmp      = ( reshapeg( c(:,i),k,k )*s )';
       impulse  = [ impulse tmp(:) ];

    end 

	impulse = impulse';
    
    % Compute cumulative sum of impulse responses for levels 
	impulse = cumsum(impulse);		

    % Construct variance decompositions 
    tmp0   = reshapeg( cumsum( impulse(1:end-1,:).^2 ),q*k,k );
    tmp1   = repmat( sum(tmp0,2),1,k );
    decomp = reshapeg( 100*tmp0 ./ tmp1 , q , k^2 );
        
    
    % Plot impulse responses 
    figure(1);
   
    subplot(5,5,1)
    plot(seqa(0,1,q+1),[impulse(1:q+1,1) zeros(q+1,1)] )
	title('Agg. supply shock');
	ylabel('Real Output');
    xlabel ('Quarter');

    subplot(5,5,2)
    plot(seqa(0,1,q+1),[impulse(1:q+1,2) zeros(q+1,1)] );
	title('Agg. demand shock');
	ylabel('Real Output');
    xlabel ('Quarter');

    subplot(5,5,3)
    plot(seqa(0,1,q+1),[impulse(1:q+1,3) zeros(q+1,1)] );     
	title('Aust. portfolio shock');
	ylabel('Real Output');
    xlabel ('Quarter');
 
    subplot(5,5,4)
	plot(seqa(0,1,q+1),[impulse(1:q+1,4) zeros(q+1,1)]);   
	title('Nominal shock');
	ylabel('Real Output');
    xlabel ('Quarter');
    
    subplot(5,5,5)
    plot(seqa(0,1,q+1),[impulse(1:q+1,5) zeros(q+1,1)]);   
	title('US portfolio shock');
	ylabel('Real Output');
    xlabel ('Quarter');

    
    %------------
    
    subplot(5,5,6)
    plot(seqa(0,1,q+1),[impulse(1:q+1,6) zeros(q+1,1)] )
	title('Agg. supply shock');
	ylabel('Interest rate');
    xlabel ('Quarter');
    
    subplot(5,5,7)
    plot(seqa(0,1,q+1),[impulse(1:q+1,7) zeros(q+1,1)] );
	title('Agg. demand shock');
	ylabel('Interest rate');
    xlabel ('Quarter');
    
    subplot(5,5,8)
    plot(seqa(0,1,q+1),[impulse(1:q+1,8) zeros(q+1,1)] );     
	title('Aust. portfolio shock');
	ylabel('Interest rate');
    xlabel ('Quarter');
    
    subplot(5,5,9)
	plot(seqa(0,1,q+1),[impulse(1:q+1,9) zeros(q+1,1)]);   
	title('Nominal shock');
	ylabel('Interest rate');
 
    subplot(5,5,10)
    plot(seqa(0,1,q+1),[impulse(1:q+1,10) zeros(q+1,1)]);   
	title('US portfolio shock');
	ylabel('Interest rate');
    xlabel ('Quarter');
    

    % ----------------

    subplot(5,5,11)
    plot(seqa(0,1,q+1),[impulse(1:q+1,11) zeros(q+1,1)] )
	title('Agg. supply shock');
	ylabel('Real Aust. equity');
    xlabel ('Quarter');
    

    subplot(5,5,12)
    plot(seqa(0,1,q+1),[impulse(1:q+1,12) zeros(q+1,1)] );
	title('Agg. demand shock');
	ylabel('Real Aust. equity');
    xlabel ('Quarter');
    

    subplot(5,5,13)
    plot(seqa(0,1,q+1),[impulse(1:q+1,13) zeros(q+1,1)] );     
	title('Aust. portfolio shock');
	ylabel('Real Aust. equity');
    xlabel ('Quarter');
    
 
    subplot(5,5,14)
	plot(seqa(0,1,q+1),[impulse(1:q+1,14) zeros(q+1,1)]);   
	title('Nominal shock');
	ylabel('Interest rate');
 
    subplot(5,5,15)
    plot(seqa(0,1,q+1),[impulse(1:q+1,15) zeros(q+1,1)]);   
	title('US portfolio shock');
	ylabel('Real Aust. equity');
    xlabel ('Quarter');
    
    % ----------------

    subplot(5,5,16)
    plot(seqa(0,1,q+1),[impulse(1:q+1,16) zeros(q+1,1)] )
	title('Agg. supply shock');
	ylabel('Price');
    xlabel ('Quarter');
    

    subplot(5,5,17)
    plot(seqa(0,1,q+1),[impulse(1:q+1,17) zeros(q+1,1)] );
	title('Agg. demand shock');
	ylabel('Price');
    xlabel ('Quarter');
    
    subplot(5,5,18)
    plot(seqa(0,1,q+1),[impulse(1:q+1,18) zeros(q+1,1)] );     
	title('Aust. portfolio shock');
	ylabel('Price');
    xlabel ('Quarter');
    
 
    subplot(5,5,19)
	plot(seqa(0,1,q+1),[impulse(1:q+1,19) zeros(q+1,1)]);   
	title('Nominal shock');
	ylabel('Price');
    xlabel ('Quarter');
    
 
    subplot(5,5,20)
    plot(seqa(0,1,q+1),[impulse(1:q+1,20) zeros(q+1,1)]);   
	title('US portfolio shock');
	ylabel('Price');
    xlabel ('Quarter');
    
    % ----------------

    subplot(5,5,21)
    plot(seqa(0,1,q+1),[impulse(1:q+1,21) zeros(q+1,1)] )
	title('Agg. supply shock');
	ylabel('Real US equity');
    xlabel ('Quarter');
    

    subplot(5,5,22)
    plot(seqa(0,1,q+1),[impulse(1:q+1,22) zeros(q+1,1)] );
	title('Agg. demand shock');
	ylabel('Real US equity');
    xlabel ('Quarter');
    
    subplot(5,5,23)
    plot(seqa(0,1,q+1),[impulse(1:q+1,23) zeros(q+1,1)] );     
	title('Aust. portfolio shock');
	ylabel('Real US equity');
    xlabel ('Quarter');
    
    subplot(5,5,24)
	plot(seqa(0,1,q+1),[impulse(1:q+1,24) zeros(q+1,1)]);   
	title('Nominal shock');
	ylabel('Real US equity');
    xlabel ('Quarter');
    
    subplot(5,5,25)
    plot(seqa(0,1,q+1),[impulse(1:q+1,25) zeros(q+1,1)]);   
	title('US portfolio shock');
	ylabel('Real US equity');
    xlabel ('Quarter');
       
    
end
%
%--------------------------- Functions ----------------------------------
% 
%--------------------------------------------------------------------------
% Log-likelihood function for restricted SVAR
%--------------------------------------------------------------------------
function lnl = neglogl( b,v,a1 )

    [ t,n ] = size( v );
    
    lf  = zeros( t,1 );
    f   =     [   b(1)   0       0       0      0       ;
        	      b(2)   b(3)    b(4)    0      0       ;
	              b(5)  -b(3)    b(6)    0      b(12)   ;
   		          b(7)   b(8)    b(9)    b(10)  0       ;
			       0     0       0       0      b(11) ] ;

	s      = a1*f;		% Long-run restrictions
	omegav = s*s';

    for i = 1:t;

        lf(i) = -0.5*n*log(2*pi) - 0.5*log(det(omegav)) - 0.5*v(i,:)*inv(omegav)*v(i,:)';    
    end
    
    lnl = -mean( lf );

end

%--------------------------------------------------------------------------
% Log-likelihood function for unrestricted SVAR
%--------------------------------------------------------------------------

function lnl = neglogu( b,v,a1 )


    [ t,n ] = size( v );
    
    lf  = zeros( t,1 );
    f   =     [   b(1)   0       0       0      0       ;
        	      b(2)   b(3)    b(4)    0      0       ;
	              b(5)   b(13)   b(6)    0      b(12)   ;
   		          b(7)   b(8)    b(9)    b(10)  0       ;
			       0     0       0       0      b(11) ] ;

	s      = a1*f;		% Long-run restrictions
	omegav = s*s';

    for i = 1:t;

        lf(i) = -0.5*n*log(2*pi) - 0.5*log(det(omegav)) - 0.5*v(i,:)*inv(omegav)*v(i,:)';    
    end
    lnl = -mean( lf );
end