%=========================================================================
%
%   Program to generate the sampling distribution of the eigenvalues of
%   a bivariate vecm with rank r=1
%
%=========================================================================
clear all
clc

format short g
RandStream.setGlobalStream( RandStream('mt19937ar','seed',1234) );
    
Tv   = [100,200,400,800];
reps = 10000;
lam  = zeros(reps,2);
lr   = zeros(reps,length(Tv));

for tc = 1:length(Tv)
    
    t = Tv(tc);

    for rep = 1:reps
            
        % Simulate dgp
        y2 = cumsum(randn(t,1));     
        y1 = y2 + randn(t,1);
        y  = [ y1  y2 ];

        r0 = trimr(y,1,0)-trimr(y,0,1);     
        r1 = trimr(y,0,1);

        % Construct sums of squares matrices     
        tmp = length(r0);
        s00 = r0'*r0/tmp;                  
        s11 = r1'*r1/tmp;
        s01 = r0'*r1/tmp;
        s10 = s01';

        % Choleski decomposition of s11         
        l = chol(s11)'; 
            
         % Compute eigenvalues          
         lam(rep,:) = eig( inv(l)*s10*inv(s00)*s01*inv(l') );
         lam(rep,:) = flipud(lam(rep,:)')';
   
         % Compute trace statistic (smallest eigenvalue is zero)  
         lr(rep,tc) = -t*log(1-lam(rep,2));       

    end

    disp('           T       Mean 1      Mean 2       StDev1      StDev2');
    disp([t mean(lam) std(lam)]);
    disp( ' ')
end


%  Compute and plot the kernel density estimator of the lr statistic 
%  and compare with chisq dof=1   
  
minx = 0.5; 
maxx = 12;
xi = seqa(minx,(maxx-minx)/200,201)';

fchi = chi2pdf(xi,1);

% Pick last column of lr corresponding to t=800 to compute kernel estimate    
fhat = ksdensity(lr(:,tc),xi);

%**********************************************************************
%***
%***     Generate graph
%***
%**********************************************************************

figure(1)
clf

plot(xi,fchi,'--k',xi,fhat,'-k');
ylabel('f(LR)');
xlabel('LR');
box off
