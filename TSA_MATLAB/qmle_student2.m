%========================================================================
%
%   Program to simulate data from a t - distribution and compute
%   scaled covariance matrices for true and misspecified models
%
%========================================================================

function qmle_student2( )

    clear all
    clc
   
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',12) )


    mu  = 10;     % Population mean      
    sig = 1;     % Population standard deviation  
    gam = 10;    % Degrees of freedom
    t   = 5000;

    % Generate data from the true model (Student t)
    v = tinv( rand(t,1),gam );           
    y = mu + sig*sqrt( (gam-2)/gam )*v;


    % ------------------ CORRECT SPECIFICATION ------------------------
    % Estimate the Student t model  
    bstart = [ mu  sig^2  gam ];


    % Optimization settings
    options = optimset( 'LargeScale','off',  ...
                        'Display','iter',     ...
                        'MaxIter',15000,     ...
                        'MaxFunEvals',10000, ...
                        'TolFun',1e-12,       ...
                        'TolX',1e-12);
    

    %[bhat,fu] = fminsearch( @(b) lnl(b,y),bstart,options); 
     [bhat,fu,~,~,~,H] = fminunc( @(b) lnl( b,y ),bstart,options);

    
    % Compute gradient and Hessian
    g = numgrad( @lnlt,bhat',y );
 
    
    iH = inv( H);
    
    format short
    
    disp('Results based on the correct specification');
    disp('------------------------------------------');
    
    disp('Covariance matrix (Hessian)');
    disp('---------------------------');
    disp( (1/t)*iH );

    
    j = g'*g/t;
    disp('Covariance matrix (OPG)');
    disp('---------------------------');
    disp( (1/t)*inv( j ) );

    
    disp('Covariance matrix (QMLE)');
    disp('---------------------------');
    disp( (1/t)*( iH*j*iH ) );

    % ------------------ INCORRECT SPECIFICATION ------------------------
    % Idea here is that the QMLE estimates should approximate the 
    % correct estimates for the mean and the variance. So compare this 
    % (2x2) matrix with the upper (2x2) submatrix of the correct model.
    
    
    % Estimate the parameters of the normal 
 
    m  = mean(y);
    s2 = mean( (y - m).^2 );
 
 
    % Compute gradients of the misspecified model (normal)
    g1 = (y - m)/s2;
    g2 = -0.5/s2 + 0.5*(y - m).^2/s2^2;
    g  = [ g1 g2 ];
    j  = g'*g/t;
 
     
    % Compute hessian of the misspecified model (normal) 
    H = zeros( 2,2 );
    H(1,1) = -1/s2;
    H(1,2) = -mean( y - m )/s2^2;
    H(2,1) = H(1,2);
    H(2,2) = 0.5/s2^2 - mean( (y - m).^2 )/s2^3;
 
    iH = -inv(H);

    disp('Results based on the incorrect specification');
    disp('------------------------------------------');
    
    disp('Covariance matrix (Hessian)');
    disp('---------------------------');
    disp( (1/t)*iH );

    disp('Covariance matrix (OPG)');
    disp('---------------------------');
    disp( (1/t)*inv( j ) );

    
    disp('Covariance matrix (QMLE)');
    disp('---------------------------');
    disp( (1/t)*( iH*j*iH ) );


end

%=======================================================================
%
%   Wrapper function      
%
%=======================================================================

function logl = lnl( b,y )

    logl = -mean( lnlt( b,y ) );

end


%=======================================================================
%
%   Log-likelihood funciton for a Student t disturbance     
%
%=======================================================================

function loglt = lnlt( b,y )


    u   = y - b(1);         
    s2  = abs( b(2) );      
    gam = abs( b(3) );                                                   

    z     = u./sqrt(s2);                                                 
    const = gamma( (gam+1)/2 ) / ( sqrt(pi*(gam-2)) * gamma( gam/2 ) );  
    loglt   = log(const) - 0.5*log(s2) - 0.5*(gam+1)*log( 1 + (z.^2)/(gam-2) );                                                                                                                                                                                        

end 

