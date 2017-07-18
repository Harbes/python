%=========================================================================
%
%   Program to compare single equation cointegration estimators
%
%=========================================================================
function coint_reg( )

    clear all
    clc
    
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) );
    
    reps = 10000;
    Tv   = [100,200,400];
    rhov = [0,0.4,0.8];
    av   = [0,0.4,0.8];

    bhat  = zeros(reps,3); 

    for ac = 1:length(av) 

        for rhoc = 1:length(rhov) 

            for tc = 1:length(Tv) 

                t = Tv(tc);

                for rep = 1:reps 

                    u      = randn(t,2); 
                    u(:,2) = rhov(rhoc)*u(:,1)+sqrt(1-rhov(rhoc)^2)*u(:,2); 
                    u(:,1) = recserar(u(:,1),u(1,1),av(ac));
        			u(:,1) = u(:,1)/std(u(:,1));

                    x = cumsum(u(:,2));
                    y = x + u(:,1);
			
                    % OLS
                    bhat(rep,1) = x\y;
                    e           = y-x*bhat(rep,1);

                    % FMOLS
                    dx = trimr(x,1,0)-trimr(x,0,1);
                    [Omega,Delta] = lrvars( [ trimr(e,1,0) dx ]);
                    ys            = trimr(y,1,0) - dx/Omega(2,2)*Omega(2,1);
                    Deltas        = Delta(1,2)-Delta(2,2)/Omega(2,2)*Omega(2,1);
                    bhat(rep,2)   = (ys'*trimr(x,1,0)-t*Deltas)/(trimr(x,1,0)'*trimr(x,1,0));


                    % VECM
                    bhat(rep,3) = vecm([y x]);

                end
                disp(['Sample size       = ',num2str(t) ])
                disp(['alpha             = ',num2str(av(ac)) ])
                disp(['rho               = ',num2str(rhov(rhoc)) ])
                
                disp('           Bias                        Std Dev           ')
                disp('   --------------------------    --------------------------')
                disp('      OLS      FMOLS      VECM      OLS      FMOLS      VECM');

                disp( [mean(bhat - 1)   std(bhat)])

            end

        end

    end

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Long-run variance
%-------------------------------------------------------------------------
function [Omega,Delta] = lrvars(z)

    [t,k] = size(z); 
    p     = floor(8*(t/100)^(2/9));
    

    J0 = z'*z; 
    J1 = 0;

    for j =1:p
        
        Gj = trimr(z,j,0)'*trimr(z,0,j);
        J0 = J0 + Gj + Gj';
        J1 = J1 + 2*j*Gj;
    end
    i = ones(k,1);
    v0 = i'*J0*i; 
    v1 = i'*J1*i;
    p = min( [ floor(1.1447*((v1/v0)^2*t)^(1/3)); p ] );

    Omega = z'*z/t; 
    Delta = z'*z/t;

    for j = 1:p
        
        Gj    = trimr(z,j,0)'*trimr(z,0,j)/t;
        Omega = Omega + (1-j/(p+1))*(Gj + Gj');
        Delta = Delta + (1-j/(p+1))*Gj;
    end
end
%-------------------------------------------------------------------------
%  Johansen VECM estimator
%-------------------------------------------------------------------------
function bhat = vecm(y)

    [t,k] = size(y); 
    pmax  = 12; 

    logL = []; 
    x    = []; 
    y0   = trimr(y,pmax,0); 

    tmp  = length(y0);
    logL = -0.5*tmp*log(det(y0'*y0/tmp));

    for j = 1:pmax 
  
        x    = [ x trimr(y,pmax-j,j) ];
        vhat = y0-x*(x\y0);
        logL = [logL ;  -0.5*length(vhat)*log(det(vhat'*vhat/length(vhat))) ];
    end
        
    nparams = k+k^2*seqa(0,1,pmax+1)';
    HQ      = -2*logL/tmp+2*nparams*log(log(tmp))/tmp;
    [~,p]   = min(HQ);
    p       = p-1;

    [~,beta,~,~,~] = johansen(y,p,1,1);
    
    bhat = -beta(2);
    
end

%-------------------------------------------------------------------------
%  Johansen procedure 
%-------------------------------------------------------------------------
function [alpha,beta,logl,maxt,tracet] = johansen(y,p,r,model)

    ty = length(y);
    dy = trimr(y,1,0)-trimr(y,0,1); 
    z0 = trimr(dy,p-1,0);
    z1 = trimr(y,p-1,1);
   
               
    z2 = [];
    
    for j =1:p-1
         
        z2 = [ z2 trimr(dy,p-1-j,j) ];
    end
    
    if model == 1;
                                      
        z1 = z1;
        z2 = z2;
        
    elseif model == 2 
            
        z1 = [ trimr(y,p-1,1)  ones(ty-p,1) ];
        z2 = z2;
        
    elseif model == 3 
            
        z1 = z1;
        z2 = [ ones(ty-p,1)  z2 ];   
     
    elseif model == 4
            
        z1 = [z1  seqa(1,1,ty-p)'];
        z2 = [ ones(ty-p,1)  z2 ];   

    elseif model == 5
            
        z1 = z1;
        z2 = [ones(ty-p,1)  seqa(1,1,ty-p)'  z2];   
    end
        
    if p == 1 && model <= 2
        
         r0 = z0; 
         r1 = z1;
         
    else
        
        r0 = z0-z2*(z2\z0); 
        r1 = z1-z2*(z2\z1);
    end
          
    [ tr,tc ]  = size( r0 );
    
    % Construct sums of squares matrices
    s00 = r0'*r0/tr;                                        
    s11 = r1'*r1/tr;
    s01 = r0'*r1/tr;
    s10 = s01';
    
    % Solve eigenvalue problem
    l       = chol(s11)';                                               
    [e,tmp] = eig( inv(l)*s10*inv(s00)*s01*inv(l') );
           
    % Sort eigenvalues and store index
    [ lam,IX ] = sort( diag(tmp),'descend' ); 

    % Sort eigenvectors on eigenvalue index and normalize
    gam = (inv(l')*e);
    gam = gam(:,IX);
   
    % Estimates
    beta  = rref(gam(:,1:r)')';                                               
    alpha = ((r1*beta)\r0)';    
    
    % Statistics
    logl   = -0.5*tc*(log(2*pi)+1) - 0.5*log(det(s00)) - 0.5*sum( log(1-lam(1:r)));
    tracet = -tr*flipud(cumsum(log(1-flipud(lam))));                     
    maxt   = -tr*log(1-lam);                  
    
end
