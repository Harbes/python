%========================================================================
%
%   Program to compute sampling properties of quasi-maximum likelihood 
%   estimator where true distribution is a negative binomial and the 
%   misspecified distribution is Poisson
%
%========================================================================
function qmle_nbar1reg( )

    clear all;
    clc;
    
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',123456) )

    Tv = [25,50,100,200,400,800,1600];

    b0    = 0.5;
    Tv    = [25,50,100,200,400,800,1600];
    nreps = 100;
    tstat = zeros(nreps,2);
    tp    = zeros(nreps,2);
    bnv   = zeros(nreps,1);
    bpv   = zeros(nreps,1);

    ops   = optimset('MaxFunEvals',20000,'MaxIter',10000);
    for k = 1:length(Tv) 
        t = Tv(k);

        for j = 1:nreps 
		
            y = ones(t+1,1);
        
            for i = 2:t+1
                % Generate negative binomial random numbers
                y(i) = nbinrnd(1+b0*y(i-1),0.5); 
            end
            x = [trimr(y,0,1) ones(t,1)]; 
            y = trimr(y,1,0);

            xxi   = inv(x'*x); 
            bhat  = xxi*x'*y;
            uhat  = y-x*bhat;
            s2    = uhat'*uhat/t; 
            seols = sqrt(s2*xxi(1,1));
            xuhat = bsxfun(@times,x,uhat);

            whitevar = xxi*(xuhat'*xuhat)*xxi; 
            seW      = sqrt(whitevar(1,1));

            % Store results
            bnv(j) = bhat(1);

            tstat(j,1) = (bhat(1)-b0)/seols;

            tstat(j,2) = (bhat(1)-b0)/seW;



%		ba = log(bhat(1)/(1-bhat(1)))|sqrt(bhat(2)); 

            ba = [log(b0/(1-b0)); 1 ]; 
            bp  = fminsearch(@(b) neglog(b,y,x),ba,ops);

            bp(1) = 1/(1+exp(-bp(1))); 
            bp(2) = bp(2)^2; 

            xdxb = bsxfun(@rdivide,x,x*bp);
            Hi   = -inv(x'*xdxb); 
            seHi = sqrt(Hi(1,1));
            g    = bsxfun(@times,xdxb,y-x*bp);

            J = g'*g;

            HiJHi = Hi*J*Hi; 
            seHJH = sqrt(HiJHi(1,1));



            bpv(j)  = bp(1);
            tp(j,1) = (bp(1)-b0)/seHi;
            tp(j,2) = (bp(1)-b0)/seHJH;


        end

        disp('   T          meanN     meanP     varN      varP      tOLS      tOLSW     tPH       tPHJH')
        tmp = [t mean([bnv bpv]) var([bnv bpv],1) mean(abs(tstat)>1.96) ...
              mean(abs(tp) > 1.96)];
        disp(tmp)  
        
          
          
    end
end
    
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
% Negative unconstrained log-likelihood
%-------------------------------------------------------------------------
function lf = neglog(b,y,x)

    b(1) = 1/(1+exp(-b(1))); 
    b(2) = b(2)^2; 
    lf   = -mean(y.*log(x*(b))-(x*(b)));
end

