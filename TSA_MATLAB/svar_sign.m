%==========================================================================
%
%      Program to estimate a svar with identification based on sign restrictions
%
%==========================================================================
function svar_sign( )

    clear all
    clc


    % GDP and CPI monthly data from 1950Q1 to 2009Q4 for U.S.
    load sign.mat         
    
    % Define variables     
    rgdp = ytdata(:,1);
    cpi  = ytdata(:,2);

    % Construct variables for use in VAR: 
    yact = [ log(rgdp)  log(cpi) ];
    y = 400*(trimr(yact,1,0) - trimr(yact,0,1));

    % Choose parameters		
    p = 2;          % Order of VAR lags		
    q = 40;         % Order of VMA lags		

    % Estimate the VAR(p) with a constant and a time trend		
    ylag = lagmatrix(y,1:p);
    ylag(any(isnan(ylag),2),:) = [];

    xmat = [ ones(length(ylag),1)  (1:1:length(ylag))'  ylag];
    ymat = trimr(y,p,0);

    bar = xmat\ymat;
    nconst = 2;
    mu = bar(1:nconst,:);

    v   = ymat - xmat*bar;
	vc = v'*v/size(v,1);

    disp('Residual variance-covariance matrix ')
    disp(vc);
    disp(' ')

    % Construct A(L) (VAR) matrices and hence the C(1) matrix 
    bar = trimr(bar,nconst,0);
	k   = size(v,2);
    a   = zeros(k^2,p);
    a1 = eye(k);

    for i = 1:p      
       		a(:,i) = vec(bar(1+k*(i-1):i*k,:));
		    a1 = a1 - reshapeg(a(:,i),k,k);
    end

    % Generate S matrix from Choleski decomposition   
    s = chol(vc)';        

    % Compute impulse responses based on Choleski 
    % Identification is based on short-run restrictions without sign restrictions    
    impulse_nosign = irf(a,s,1,k,q,p);     

    disp('IRF based on short-run restrictions without sign restrictions')
    disp(impulse_nosign);
    disp(' ');

    % Compute impulse responses rotated by the matrix Q'    
    th   = 0.2*pi;
    qmat =[ cos(th) , -sin(th) ;
           sin(th) ,  cos(th) ];

    impulse_rotate = irf(a,s*qmat',1,k,q,p);     

    disp('IRF based on based on Choleski which is rotated by the matrix Q')
    disp(impulse_rotate);
    disp(' ');



    % Select IRFs that satisfy sign restrictions      
    nsim = 10000;         % Number of draws to compute IRFs    
    impulse_select = zeros(q+1,1);

    nsuccess = 0.0;
    for iter = 1:nsim
        th = pi*rand(1,1);
        qmat = [cos(th) , -sin(th);
                sin(th) ,  cos(th)   ] ;

        impulse = irf(a,s*qmat',1,k,q,p);      
    
        % Choose impulses that satisfy the sign restrictions for all values
        if min(min(impulse(:,[1 2 4]))') >= 0.0 && min(min(impulse(:,3))') < 0.0          

         impulse_select = [impulse_select , impulse];

         nsuccess = nsuccess + 1;

        end
    end

	disp(['Number of simulations          = ', num2str(nsim) ]);
    disp(['Number of successful draws (%) = ', num2str(100*nsuccess/nsim) ]);

    % Find the model that corresponds to the median impulses         
    impulse_select = trimr(impulse_select',1,0)';  
    
     % Choose contemporaneous impulses = nsuccess x cols(y)^2
    impulse_contemp = reshape(impulse_select(1,:),size(impulse_select,2)/size(y,2)^2,size(y,2)^2);         

    zimpulse_contemp = (impulse_contemp - repmat(median(impulse_contemp),size(impulse_contemp,1),1))./repmat(std(impulse_contemp),size(impulse_contemp,1),1);        %     Compute deviations  

    total = sum(zimpulse_contemp')';                                                                %     Compute total deviations for each model     

    [~,kmin] = min(abs(total));                                                                     %     Choose the set of impulses which yields the minimum absolute total deviation for each model  

    impulse_med = impulse_select(:,((kmin-1)*size(y,2)^2+1):(kmin*size(y,2)^2));                        %     Overall median impulse      

    % Plot IRFs
    figure(1)	
    subplot(2,2,1);
	plot(seqa(0,1,q+1),impulse_med(1:q+1,1)); 
    title('Supply shock');
	ylabel('Output');
    xlabel('Quarter');

	subplot(2,2,2);
	plot(seqa(0,1,q+1),impulse_med(1:q+1,2));   
	title('Demand shock');
	ylabel('Output');

	subplot(2,2,3);
	plot(seqa(0,1,q+1),impulse_med(1:q+1,3));   
	title('Supply shock');
	ylabel('Price');

	subplot(2,2,4);
	plot(seqa(0,1,q+1),impulse_med(1:q+1,4));   
	title('Demand shock');
	ylabel('Price');

end
%--------------------------- Functions ----------------------------------
% 
%--------------------------------------------------------------------------
%   Compute impulse responses
%   itype=0 for non-accumulated and itype=1 for accumulated
%--------------------------------------------------------------------------
function impulse = irf(phi,s,itype,k,q,p)

    % Construct SI(L) matrices (VMA)     
    si = vec(eye(k));

    for i = 1:q

        sum = 0;
        for j = 1:min([p i]);
            sum = sum + reshapeg(phi(:,j),k,k)*reshapeg(si(:,i-j+1),k,k);
        end
        si  = [si , vec(sum')];
    end

    % Construct orthogonal impulse response functions        
    impulse = vec(s');

    for i = 2:q+1
       impulse  = [impulse , vec( (reshapeg(si(:,i),k,k)*s )' )];
    end
	impulse = impulse';

    if itype
        impulse = cumsum(impulse);  
    end                 %		Compute cumulative sum of impulse responses for levels

end


