%=========================================================================
%
%   Program to estimate impulse responses for a VECM and a VAR
%   of money demand
%
%=========================================================================
function coint_impulse( )

    clear all
    clc
    
    % Load the data
    load moneydemand.mat

    lrm = log(m2./cpi);
    lry = log(gdp./cpi);
    ym  = [ lrm  lry  tbill/100  fedfunds/100 ];

% // Variable names and horizon for impulse responses
% vars = "Money" $| "Income" $| "T-Bill" $| "Fed Funds";
% horizon = 40;

    t = 188;
    y = ym(1:t,:);

    p     = 4;      % Lags in VAR                    
    r     = 1;      % Number of cointegrating equations   
    model = 3;
    
    % Impose unit coefficient on y
    H = [ 1 0 0;
         -1 0 0;
          0 1 0;
          0 0 1];

    [ lam,~,beta,~,~,~ ] = johansenH(y,p,r,model,H);
    betahat = H*beta; 
  
    % Estimate remaining VECM coefficients, model 3
    dy = trimr(y,1,0)-trimr(y,0,1); 
    x = [ones(t-p,1) trimr(y,p-1,1)*betahat ];

    for j = 1:p-1 

        x = [x trimr(dy,p-1-j,j)];
    end
    bhat = x\trimr(dy,p-1,0);
    vhat = trimr(dy,p-1,0)-x*bhat;
    
    ghat = trimr(bhat,r+1,0);
    
    c    = size(y,2);
    gam1 = eye(c);
    
    for j = 1:p-1
        
        gam1 = gam1 - ghat((j-1)*c+1:j*c,:)';

    end
   
    % Put VECM in levels VAR form  
    alphahat = bhat(2:r+1,:)';
    ghat     = [trimr(bhat,r+1,0)' zeros(c,c) ];
    A        = eye(c)+alphahat*betahat'+ghat(:,1:c);
    
    for j = 2:p

        tmp1 = (j-1)*c;
        tmp2 = (j-2)*c;
        A = [ A ghat(:,tmp1+1:j*c)-ghat(:,tmp2+1:tmp1) ];
        
    end
 
    vecm_impulse = impulse(A,chol(vhat'*vhat/length(vhat))',40,c);
    
    %  Estimate levels VAR for comparison with VECM
    x = ones(t-p,1);

    for j = 1:p

        x = [x  trimr(y,p-j,j)];

    end
    bhat = x\trimr(y,p,0);
    vhat = trimr(y,p,0)-x*bhat;

    var_impulse = impulse(trimr(bhat,1,0)',chol(vhat'*vhat/length(vhat))',40,c);
       
    % VECM impulse responses to interpret coefficient on income
    k = [0;0;0;0];
    k = [k zeros(4,3)]; 

    k(2,2) = 1;
    k(1,2) = 1;

    vecm1_impulse = impulse(A,gam1*k,40,c);
    
    hh = 0:1:40;
    display('   Illustrating the unit long run elasticity of money wrt income' )
    display('----------------------------------------------------------------')
    display('    Time      Money     Income    Tbill    FedFunds   ' )
    disp([hh' vecm1_impulse(:,5:8) ]);

    
    figure(1)
    subplot(3,1,1)
    plot(hh,var_impulse(:,5:8))
    title('Impulses associated with an income shock: VAR in levels')

    subplot(3,1,2)
    plot(hh,vecm_impulse(:,5:8))
    title('Impulses associated with an income shock: VECM')
    
    subplot(3,1,3)
    plot(hh,vecm1_impulse(:,5:8))
    title('Impulses associated with an income shock: VECM with unit restriction')
end

%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Johansen procedure 
%-------------------------------------------------------------------------
function [lam,alpha,beta,logl,tracet,maxt] = johansenH(y,p,r,model,H)

    ty = length(y);
    dy = trimr(y,1,0)-trimr(y,0,1); 
    z0 = trimr(dy,p-1,0);
    z1 = trimr(y,p-1,1)*H;
              
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
            
        z1 = [z1  seqa(1,1,ty-p)];
        z2 = [ ones(ty-p,1)  z2 ];   

    elseif model == 5
            
        z1 = z1;
        z2 = [ones(ty-p,1)  seqa(1,1,ty-p)  z2];   
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
    logl   = 0.5*(-tr*tc*(log(2*pi)+1) - tr*log(det(s00)) - tr*sum( log(1-lam(1:r))));
    tracet = -tr*flipud(cumsum(log(1-flipud(lam))));                     
    maxt   = -tr*log(1-lam);                  
    
end

%-------------------------------------------------------------------------
%  Impulse responses from VAR(p) 
%  y(t) = det + A1 y(t-1) + ... + A2 y(t-p) + v(t)
%  A is coefficient matrix:  A = ( A1 A2 ... Ap 
%  C contains shocks (eg chol of variance matrix of v(t)
%  n is dimension of y
%-------------------------------------------------------------------------
function imp = impulse(A,C,h,n)

    [ rA,cA ] = size(A);
    [ rC,cC ] = size(C);

    p = cA/rA;

    %  VAR(1) companion form
    A = [ A ;
          [ eye((p-1)*n) zeros((p-1)*n,n) ] ];

    C = [ [ C  zeros(rC,(p-1)*cC) ] ; zeros((p-1)*cC,p*cC) ];

    % Impulse responses
    
    AC = C;

    tmp = AC(1:n,1:n);
    imp = tmp(:)';
 
    for j = 1:h
        
        AC  = A*AC;
        tmp = AC(1:n,1:n);
        imp = [ imp; tmp(:)'];
        
    end


end

%     % Compute orthonormal complement of alpha
%     alpha0 = null(alpha');      
%     %disp(['Orthogonal complement matrix alpha0 = ', num2str(alpha0') ]);
% 
%     % Computer constant term in cointegrating regression
%     beta0 = (alpha'*alpha)\(alpha'*mu0);
% 

