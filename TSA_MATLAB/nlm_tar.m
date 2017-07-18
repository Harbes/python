%=========================================================================
%
%   Nonparametric density estimator of a threshold autoregressive model
%
%=========================================================================

function nlm_tar( )

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',3) );

    theta = 0.8; 
    k = 1;                  
    t = 5000;
    y = -3:0.1:3;


    % Simulate TAR(1) model
    yt = theta*zeros(t,1);

    for i = 2:t

        yt(i) = theta*abs(yt(i-1)) + sqrt(1 - theta^2)*randn(1,1);

    end


    % Kernel regression of y_t on y_t-k 
    
    fx  = zeros(length(y),1);
    fxy = zeros(length(y),1);

    h = 1.06*std(yt)*t^(-1/5);

    
    for i = 1:length(y);
        z      = ((y(i) - trimr(yt,0,k))./h);    
        fx(i)  = mean( normpdf(z)./h );
        fxy(i) = mean( normpdf(z).*trimr(yt,k,0)./h );
    end

    m = fxy ./ fx;

    % Compute true conditional mean
    m_true = theta*abs(y);


    % Compute linear conditional mean 
    m_linear = y*theta.^k;

    %**********************************************************************
    %***
    %***     Generate graph
    %***
    %**********************************************************************

    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;

    %--------------------------------------------------------%
    plot(y,m,'--k','LineWidth',1.0);
    hold on;
    plot(y,m_true,'-k','LineWidth',0.5);
    plot(y,m_linear,':k','LineWidth',1.00);
    hold off;
    ylabel('$m(y_{t-1})$');
    xlabel('$y_{t-1}$');
    % legend('Estimated','True','Location','NorthWest')
    % legend boxoff;
    box off;
    hold off;
    
    %laprint(1,'nprtar','options','factory');

    
end



