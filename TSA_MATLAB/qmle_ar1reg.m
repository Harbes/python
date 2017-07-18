%=========================================================================
%
%   Sampling properties of the White estimator 
%
%=========================================================================

clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',123456) )

b0    = 0.5;
Tv    = [25,50,100,200,400,800,1600];
nreps = 10000;
tstat = zeros(nreps,2);
bhatv = zeros(nreps,1);

for j = 1:length(Tv)
    t = Tv(j);
    for k = 1:nreps
		v = randn(t+1,1);
		y = v;
        for i = 2:t+1
			y(i) = b0*y(i-1) + sqrt(1+0.5*y(i-1).^2)*v(i);
        end
		x     = trimr(y,0,1); 
        y     = trimr(y,1,0);
		xxi   = inv(x'*x); 
		bhat  = xxi*x'*y;
		uhat  = y-x*bhat;
		s2    = uhat'*uhat/t; 
        seols = sqrt(s2*xxi(1,1));
		xuhat = x.*uhat;

        % White's estimate of variance
		whitevar = xxi*(xuhat'*xuhat)*xxi; 
        seW      = sqrt(whitevar(1,1));

        % Save results
		bhatv(k)   = bhat(1);
		tstat(k,1) = (bhat(1)-b0)/seols;
		tstat(k,2) = (bhat(1)-b0)/seW;


    end

    disp(['Results for T = ', num2str(t) ]);
    disp(' Mean(bhat)  Var(bhat)   OLS    White')
	disp( [mean(bhatv) std(bhatv)^2 mean(abs(tstat) > 1.96) ]);

end