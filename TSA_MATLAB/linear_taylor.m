%=========================================================================
%
%   Program to estimate parameters of a Taylor Rule
%
%=========================================================================

function linear_taylor( )

    clear all;
    clc;

    load taylor;

    % Choose data 1987:Q1 to 1999:Q4
    infl = taylor.Data(104:end,1);  
    ygap = taylor.Data(104:end,2);
    ffr  = taylor.Data(104:end,3);

    t = length(ffr);
        
    tmp = [ t,         sum(infl),       sum(ygap),       ;
            sum(infl), sum(infl.*infl), sum(infl.*ygap)  ;
            sum(ygap), sum(infl.*ygap), sum(ygap.*ygap) ]; 
    
    disp( tmp );
    
    tmp1 = [    sum(ffr)       ;
                sum(ffr.*infl) ;
                sum(ffr.*ygap) ];
    
    disp( tmp1 )
        
    % OLS estimates long-hand
    betahat1 = inv(tmp)*tmp1;
    
    % OLS estimates
    x = [ones(t,1)  infl  ygap];
    y = ffr;
    betahat  = x\y;
    
    disp( 'ML/OLS estimates - both methods' );
    disp( [betahat1 betahat ] );
    
    e       = y - x*betahat;              
    sig2hat = e'*e/t;
    disp( 'ML/OLS estimate of variance' );
    disp( sig2hat );
   
    vcov = sig2hat*inv(x'*x);
    disp( 'Covariance matrix of parameters' )
    disp( vcov );
     
    % Wald test of restrictions b(2)=1.5 and b(3)=0.5
    R = [0 1 0 ;
         0 0 1];
    Q = [ 1.5  ; 
          0.5 ];

    W = (R*betahat - Q)'*inv(R*vcov*R')*(R*betahat - Q);
    p = 1 - chi2cdf(W,2);

    disp('Wald test results')
    disp(['Wald statistic    = ' num2str(W) ]);
    disp(['p value           = ' num2str(p) ]);

end
