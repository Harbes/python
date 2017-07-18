%=========================================================================
%
%   Program to approximate asymptotic critical values for the trace test
%
%=========================================================================
clear all
clc

RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234) );

model = 5;                  % Choose model (1 - 5)                  
t     = 1000;                   
nreps = 100000;             
tr    = zeros(nreps,1);        

pv = [ 0.90, 0.95,  0.99];    

q = [];


for nr = 1:6                %  Up to 6 common trends   

  for j = 1:nreps 

    db = randn(t,nr); 
    f  = cumsum(db);

    if model==2
        
        f = [f  ones(t,1)];

    elseif model==3 
        
        f(:,1) = seqa(1,1,t)'; 
        f      = bsxfun(@minus,f,mean(f));

    elseif model==4
        
        f = [f  seqa(1,1,t)'];  
        f = bsxfun(@minus,f,mean(f));

    elseif model==5
            
        f(:,1) = (seqa(1,1,t).^2)'; 
        x      = [ones(t,1) seqa(1,1,t)']; 
        f      = f - x*(x\f);

    end

    %  Note that db is computed as a forward difference relative to f     
	m1 = trimr(f,0,1)'*trimr(db,1,0)/t;      
	m2 = trimr(f,0,1)'*trimr(f,0,1)/(t^2);

	tr(j) = sum( diag( m1'*inv(m2)*m1 ) );    

  end

    q = [q  quantile(tr,pv)'];

end
format short
disp(['Number of common trends   = ',num2str(model) ])
disp('-----------------------------------------------')
disp('      n-r       1         2         3         4         5         6')   
disp( [pv' q] )
 

