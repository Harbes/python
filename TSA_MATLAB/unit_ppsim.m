%=========================================================================
%
%   Monte Carlo experiment on the PP coefficient tes
%
%=========================================================================
function unit_ppmc( )

    clear all
    clc

    reps = 10000;
    Tv   = [100,200];
    arv  = [0,0.3,0.6,0.9];
    mav  = [-0.8,-0.4,0,0.4,0.8];

    rej = zeros(length(arv),length(mav));
    phitilde = zeros(reps,1);

    disp( '     T          arv      mav       Rej Freq ')
    disp('---------------------------------------------')
    for i = 1:length(Tv)
                
        for j = 1:length(arv)
            
            for m = 1:length(mav)
                
                RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234) )
                
                for k = 1:reps
                    
                    v = randn(Tv(i)+1,1);
                    y = cumsum(recserar(trimr(v,1,0)+mav(m)*trimr(v,0,1),v(1),arv(j)));

                    phihat = trimr(y,0,1)\trimr(y,1,0);
                    vhat   = trimr(y,1,0)-phihat*trimr(y,0,1);
                    sig2   = vhat'*vhat/length(vhat);
                    om2    = scQS(vhat);
                    
                    phitilde(k) = phihat - 0.5*(om2-sig2)/mean(trimr(y,0,1).^2);
                end
                rej(j,m) = mean(length(y)*(phitilde-1) < -8.08);              
                disp( [Tv(i) arv(j) mav(m) rej(j,m)])
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
function om = scQS(u)

    a = trimr(u,0,1)\trimr(u,1,0);
	
    if abs(a)>0.97
        a = 0.97*a/abs(a);
    end
    
    up = trimr(u,1,0)-trimr(u,0,1)*a;
    k  = length(up);
	g  = zeros(k,1);
	
    for j = 0:k-1
		
        g(j+1) = trimr(up,j,0)'*trimr(up,0,j)/k;
		
    end

	n  = round(3*(k/100)^(2/25));
	s0 = g(1)+2*sum(g(2:n+1));
	s2 = 2*(seqa(1,1,n).^2)*g(2:n+1);

	ST  = 1.3221*(abs(s2/s0).^0.4)*k^0.2;
	x   = seqa(1,1,k-1)/ST;
	kQS = (25/(12*pi^2))./(x.^2).*(sin(6*pi*x/5)./(6*pi*x/5)-cos(6*pi*x/5));

	om = (g(1)+2*kQS*g(2:k))/(1-a)^2;
    
end