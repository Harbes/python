%=========================================================================
%
%   Program to compare the performance of alternative tests of cointegration
%
%=========================================================================

clear all
clc

RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234) );

reps = 10000;
t    = 100;                                 
n    = 2;        
r    = seqa(0,1,n+1)';                                                     
m    = n*r + (n-r).*r+ n*(n+1)/2;  % Number of free parameters                               
trcv = [12.285; 4.173];            % 5% cvs for common trends of n-r = 2,1  

% Parameters
phi1 = 0.8;                             
phi2 = 0.8;
rho  = 0.0;
phi  = diag( [phi1 phi2] );

omegav      = [ 1.0  rho  ;
                rho  1.0 ] ;
chol_omegav = chol(omegav);

% Initialise arrays
rtr = zeros(reps,1); 
rbic = zeros(reps,1); 
raic = zeros(reps,1); 
rhic = zeros(reps,1); 

for i = 1:reps 

	u  = randn(t,n)*chol_omegav;                                        
	y  = recserar(u,zeros(1,n),diag(phi)');
	
    % Construct variables for VECM
    dy = trimr(y,1,0) - trimr(y,0,1);   
	y1 = trimr(y,0,1);
    r0 = dy; 
    r1 = y1;                   
    
    % Perform eigen decomposition      
    tmp = length(r1);
	s11=r1'*r1/tmp; 
    s10=r1'*r0/tmp; 
    s00=r0'*r0/tmp;

    l = chol(s11)';                        

    lam  = eig( inv(l)*s10*inv(s00)*s10'*inv(l') );
	lam  = flipud(lam);

    [k,c] = size( r0 );
	lnl = -c*k*0.5*(log(2*pi) + 1)  - 0.5*k*log(det(s00)) - 0.5*k*([0;cumsum(log(1 - lam))] ); 


    %  Compute test statistics    
	tr      = -2*( trimr(lnl,0,1) - lnl(n+1) );     

    % LR Trace
    [~,ind] = max([tr < trcv ; 1]);
	rtr(i)  = r(ind);
    
    % AIC
    [~,ind] = min(-2*lnl + 2*m);
	raic(i) = r(ind);
    
    % BIC
    [~,ind] = min(-2*lnl + log(t)*m);
	rbic(i) = r(ind);
    
    % HIC
    [~,ind] = min(-2*lnl + 2*log(log(t))*m);
	rhic(i) = r(ind);

end

% Collect results
trace = mean(bsxfun(@eq,rtr,seqa(0,1,n+1)));
aic   = mean(bsxfun(@eq,raic,seqa(0,1,n+1)));
bic   = mean(bsxfun(@eq,rbic,seqa(0,1,n+1)));
hic   = mean(bsxfun(@eq,rhic,seqa(0,1,n+1)));

format short;

disp('     Rank    LR(trace)   BIC       AIC      HIC')
disp([r trace' aic' bic' hic'])

 

