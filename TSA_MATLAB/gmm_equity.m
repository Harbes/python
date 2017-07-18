%=========================================================================
%
%   Decompose equity returns based on GMM estimates
%   Data are daily, July 29th 2004 to March 3rd 2009, SP500,FTSE100,EURO50
%
%=========================================================================
function gmm_equity( )

    clear all
    clc

    % Read in data - 1199x3 array named 'p'
    load equity_decomposition.mat
    
    % Compute centered equity returns
    rt = 100*(log(trimr(p,1,0)) - log(trimr(p,0,1)));                          
    et = bsxfun(@minus,rt,mean(rt));  
    t  = length(et);

    % Generate squares and cross products of columns of et
    [rows,cols] = size(et);
    e2          = zeros(rows,cols*cols);
    k=1;
    for j = 1:cols   
        for i = 1:cols               
            e2(:,k) = et(:,j).*et(:,i);
            k=k+1;	
       end
    end
    
    % Drop repeated columns
    e2 = e2(:,vech( reshapeg(1:1:cols*cols,cols,cols) ) );

    % Estimate model
    ops   = optimset('LargeScale','off','Display','off');
    start = rand(6,1) ;

    [theta,qu,~,~,~,hess] = fminunc(@(b) gmmcrit(b,e2),start,ops);

    vc = inv(hess)/t;
    disp(' ');
    disp(['Value of objective function = ', num2str(qu) ]);
	disp(' ');
    disp('    Params    Std Errors '  )
    disp([theta   sqrt(diag(vc)) ] )

    disp(['J statistic             = ', num2str(2*t*qu) ]);
    nu = size(e2,2)-length(theta);
    if  nu > 0.0
        disp(['p-value             = ', num2str(1-chi2cdf(2*t*qu,nu)) ]);
    end

    disp(' ')
    disp('Covariance matrix of data ')
    disp( et'*et/t );

    a = [ theta(1:3)   diag(theta(4:6)) ] ;
    disp(' ')
    disp('Covariance matrix based on decomposition ')
    disp( a*a' )

        
    total = theta(1:3).^2 + theta(4:6).^2;
    disp(' ')
    disp('   Common   Idiosyncratic ') 
    disp( 100*( [ (theta(1:3).^2)./total (theta(4:6).^2)./total   ] ) );

end

%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
% GMM unconstrained objective function  
%-------------------------------------------------------------------------
function q = gmmcrit(b,e2)
	
    lam  = b(1:3);
    phi  = b(4:6);
   
    a        = [ lam diag(phi)] ;
    
    m = bsxfun( @minus,e2,vech(a*a')' );
    w = m'*m/length(m);
    q = 0.5*mean(m)*inv(w)*mean(m)';
end

