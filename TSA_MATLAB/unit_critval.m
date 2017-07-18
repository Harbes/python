%=========================================================================
%
%   Simulate critical values for DF, M and P tests 
%
%=========================================================================
function unit_critval( ) 

    clear all
    clc
    
    RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234) )

    reps = 100000;
    
    Tv = [ 25,50,100,250,500,1000 ];
    pv = [ 0.01,0.05, 0.10 ];

    adf1     = zeros(length(Tv),length(pv)); 
    adf2     = adf1;
    adf1cols = adf1; 
    adf2cols = adf1; 
    adf1cgls = adf1; 
    adf2cgls = adf1;
    adf1tols = adf1; 
    adf2tols = adf1; 
    adf1tgls = adf1; 
    adf2tgls = adf1;
    
    MSB     = adf1; 
    MSBcols = adf1; 
    MSBcgls = adf1; 
    MSBtols = adf1; 
    MSBtgls = adf1;
    MZtcols = adf1; 
    MZtcgls = adf1; 
    MZttols = adf1; 
    MZttgls = adf1;
    
    Pc = adf1; 
    Pt = adf1;
    
    adf = zeros(reps,2); 
    M   = zeros(reps,3);
    
    adfcols = zeros(reps,2); 
    Mcols   = zeros(reps,3);
    adfcgls = zeros(reps,2); 
    Mcgls   = zeros(reps,3);
    adftols = zeros(reps,2); 
    Mtols   = zeros(reps,3);
    adftgls = zeros(reps,2); 
    Mtgls   = zeros(reps,3);
    Pcgls   = zeros(reps,1); 
    Ptgls   = zeros(reps,1);

    for Tc = 1:length(Tv)
  
        t = Tv(Tc);
 
        for rep = 1:reps

            y     = cumsum(randn(t,1));
            ycols = glsdetrend(y,ones(t,1),-t);									 
            ytols = glsdetrend(y,[ones(t,1) seqa(1,1,t)'],-t);			 
            ycgls = glsdetrend(y,ones(t,1),-7);									 
            ytgls = glsdetrend(y,[ones(t,1) seqa(1,1,t)'],-13.5);		

            adf(rep,:)     = adftests(y,0)'; 
            M(rep,:)       = mtests(y,0)';
            adfcols(rep,:) = adftests(ycols,0)'; 
            Mcols(rep,:)   = mtests(ycols,0)';
            adfcgls(rep,:) = adftests(ycgls,0)'; 
            Mcgls(rep,:)   = mtests(ycgls,0)'; 
            adftols(rep,:) = adftests(ytols,0)'; 
            Mtols(rep,:)   = mtests(ytols,0)';
            adftgls(rep,:) = adftests(ytgls,0)'; 
            Mtgls(rep,:)   = mtests(ytgls,0)';
            Pcgls(rep)     = Ptest(y,ones(t,1),-7,0);
            Ptgls(rep)     = Ptest(y,[ones(t,1) seqa(1,1,t)'],-13.5,0);

        end

        adf1(Tc,:) = quantile(adf(:,1),pv)'; 
        adf2(Tc,:) = quantile(adf(:,2),pv)';
        adf1cols(Tc,:) = quantile(adfcols(:,1),pv)'; 
        adf2cols(Tc,:) = quantile(adfcols(:,2),pv)';
        adf1cgls(Tc,:) = quantile(adfcgls(:,1),pv)'; 
        adf2cgls(Tc,:) = quantile(adfcgls(:,2),pv)';
        adf1tols(Tc,:) = quantile(adftols(:,1),pv)'; 
        adf2tols(Tc,:) = quantile(adftols(:,2),pv)';
        adf1tgls(Tc,:) = quantile(adftgls(:,1),pv)'; 
        adf2tgls(Tc,:) = quantile(adftgls(:,2),pv)';
        MSB(Tc,:)      = quantile(M(:,2),pv)';
        MSBcols(Tc,:)  = quantile(Mcols(:,2),pv)';  
        MSBcgls(Tc,:)  = quantile(Mcgls(:,2),pv)';
        MSBtols(Tc,:)  = quantile(Mtols(:,2),pv)';  
        MSBtgls(Tc,:)  = quantile(Mtgls(:,2),pv)';
        MZtcols(Tc,:)  = quantile(Mcols(:,3),pv)';  
        MZtcgls(Tc,:)  = quantile(Mcgls(:,3),pv)';
        MZttols(Tc,:)  = quantile(Mtols(:,3),pv)';  
        MZttgls(Tc,:)  = quantile(Mtgls(:,3),pv)';
        Pc(Tc,:)       = quantile(Pcgls,pv)';   
        Pt(Tc,:)       = quantile(Ptgls,pv)';

    end  

    format short g
    
    disp('Percentiles of Distributions of Unit Root Test Statistics')
    disp('---------------------------------------------------------');

    disp('DF a test')
    disp('---------------------------------------------------------');
    disp('No de-trending');
    disp('            T         1%          5%           10% ')
    disp( [Tv' adf1] );
    disp('OLS, constant');
    disp('            T         1%          5%           10% ')
    disp( [Tv' adf1cols] );
    disp('OLS, constant and linear trend');
    disp('            T         1%          5%           10% ')
    disp( [Tv' adf1tols] );
    disp('GLS, constant (ays. equiv to no de-trending)');
    disp('            T         1%          5%           10% ')
    disp( [Tv' adf1cgls] );
    disp('GLS, constant and linear trend');
    disp('            T         1%          5%           10% ')
    disp( [Tv' adf1tgls] );

    disp('DF t test')
    disp('---------------------------------------------------------');
    disp('No de-trending');
    disp('            T         1%          5%           10% ')
    disp( [Tv' adf2] );
    disp('OLS, constant');
    disp('            T         1%          5%           10% ')
    disp( [Tv' adf2cols] );
    disp('OLS, constant and linear trend');
    disp('            T         1%          5%           10% ')
    disp( [Tv' adf2tols] );
    disp('GLS, constant (ays. equiv to no de-trending)');
    disp('            T         1%          5%           10% ')
    disp( [Tv' adf2cgls] );
    disp('GLS, constant and linear trend');
    disp('            T         1%          5%           10% ')
    disp( [Tv' adf2tgls] );
    
    
    
    disp('MZ t test')
    disp('---------------------------------------------------------');
    disp('OLS, constant');
    disp('            T         1%          5%           10% ')
    disp( [Tv' MZtcols] );
    disp('OLS, constant and linear trend');
    disp('            T         1%          5%           10% ')
    disp( [Tv' MZttols] );

    
    disp('MSB test')
    disp('---------------------------------------------------------');
    disp('No de-trending');
    disp('            T         1%          5%           10% ')
    disp( [Tv' MSB] );
    disp('OLS, constant');
    disp('            T         1%          5%           10% ')
    disp( [Tv' MSBcols] );
    disp('OLS, constant and linear trend');
    disp('            T         1%          5%           10% ')
    disp( [Tv' MSBtols] );
    disp('GLS, constant (ays. equiv to no de-trending)');
    disp('            T         1%          5%           10% ')
    disp( [Tv' MSBcgls] );
    disp('GLS, constant and linear trend');
    disp('            T         1%          5%           10% ')
    disp( [Tv' MSBtgls] );


    disp('Point optimal')
    disp('---------------------------------------------------------');
    disp('GLS, constant');
    disp('            T         1%          5%           10% ')
    disp( [Tv' Pc] );
    disp('GLS, constant and linear trend');
    disp('            T         1%          5%           10% ')
    disp( [Tv' Pt] );
 
end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Detrending function: 
%       cbar = -7 constant; 
%       cbar = -13.5 linear trend  
%       cbar = -T for OLS detrending
%-------------------------------------------------------------------------
function [ u,b ] = glsdetrend( y,x,cbar )

    t = length( y );
    yc = [ y(1)  ; (trimr(y,1,0)-(1+cbar/t)*trimr(y,0,1)) ];
    xc = [x(1,:) ; (trimr(x,1,0)-(1+cbar/t)*trimr(x,0,1)) ];
    
    b  = xc\yc;
    u  = y - x*b;
 
end
%-------------------------------------------------------------------------
%  ADF coefficient and t tests: u must be already de-trended
%-------------------------------------------------------------------------
function adfstats= adftests(u,k)


    du = trimr(u,1,0) - trimr(u,0,1);
    x  = trimr(u,0,1);

    % Set up lags of du
    if k>0;
    
        ldu = lagmatrix(du,1:1:k);
        x   = trimr( [x ldu],k,0 );
    end
    
    xxi = inv(x'*x); 
    b   = xxi*x'*trimr(du,k,0); 
    e   = trimr(du,k,0)-x*b; 
    s2  = e'*e/length(e);

    adfstats = [length(u)*b(1) ; b(1)/sqrt(s2*xxi(1,1))];

end
%-------------------------------------------------------------------------
%  M tests
%-------------------------------------------------------------------------
function tests = mtests(u,k)

    s2  = ar1rvar(u,k);
    n   = length(u);
    tmp = sum(u(1:n-1).^2);
    
    %disp(['Sum of u(1) to u(t-1) = ',num2str(tmp) ])
    %disp(['Last value: u(t)      = ',num2str(u(n)) ])
    %disp(' ') 

    u2  = tmp/n^2;     
    mza = (u(n)^2/n-s2)/(2*u2);     
    msb = sqrt(u2/s2);
    mzt = msb*mza;
    
    tests = [mza msb mzt];
 
end
%-------------------------------------------------------------------------
%  Autoregressive long run variance estimator
%-------------------------------------------------------------------------
function s2 = ar1rvar(u,k)

    du = trimr(u,1,0)-trimr(u,0,1); 
    x  = trimr(u,0,1);

     if k > 0 
        x  = [ x lagmatrix(du,seqa(1,1,k)) ];
        x(any(isnan(x),2),:) = [];
     end


    b = x\trimr(du,k,0); 
    e = trimr(du,k,0) - x*b; 
    s2 = e'*e/length(e);

    if k>0 
        s2 = s2/(1-sum(trimr(b,1,0)))^2; 
    end

end

%-------------------------------------------------------------------------
%  P-test
%-------------------------------------------------------------------------

function pt = Ptest(y,x,cbar,k)

    n  = length(y);
    uc = glsdetrend(y,x,cbar); 
    u0 = glsdetrend(y,x,0);
    s2 = ar1rvar(uc,k);

    uc = [ uc(1) ; trimr(uc,1,0)-(1+cbar/n)*trimr(uc,0,1) ];
    u0 = [ u0(1) ; trimr(u0,1,0)-trimr(u0,0,1) ];
    
    pt = (uc'*uc-(1+cbar/n)*u0'*u0)/s2;

end
