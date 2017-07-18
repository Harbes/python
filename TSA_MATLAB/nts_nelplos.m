%=========================================================================
%
%   Plot and compute ACF of the Nelson-Plosser data (1860 to 1970).
%
%=========================================================================          
function nts_nelplos( )

    clear all;
    clc;

    RandStream.setDefaultStream( RandStream('mt19937ar','seed',5) );

    % Load Nelson-Plosser data set
    load nelson_plosser.mat

    % Variable 'data' constains the following variables
    % Date 1860 - 1970
    % Real GNP, 1909 to 1970  
    % Nominal GNP, 1909 to 1970
    % Real per capita GNP, 1909 to 1970
    % Industrial production, 1860 to 1970
    % Employment, 1890 to 1970
    % Unemployment rate, 1890 to 1970
    % GNP deflator, 1889 to 1970
    % Consumer prices index, 1860 to 1970
    % Wages, 1900 to 1970
    % Real wages, 1900 to 1970 
    % Money stock, 1889 to 1970
    % Velocity, 1869 to 1970
    % Bond yield, 1900 to 1970
    % SP500, 1871 to 1970

    % Take logs except for bond yield
    rgnp   = log(data(50:111,2));
    gnp    = log(data(50:111,3));
    pcrgnp = log(data(50:111,4));
    ip     = log(data(1:111,5));
    emp    = log(data(31:111,6));
    un     = log(data(31:111,7));
    prgnp  = log(data(30:111,8));
    cpi    = log(data(1:111,9));
    wg     = log(data(41:111,10));
    rwg    = log(data(41:111,11));
    m      = log(data(30:111,12));
    vel    = log(data(10:111,13));
    bnd    = data(41:111,14);
    sp500  = log(data(12:111,15));

    % Generate deterministic and stochastic variables 
    t        = 111;
    y_dtrend = 0.1 + 0.2*(1:1:t)' + randn(t,1);      
    y_strend = recserar(0.3 + randn(t,1),0.0,1.0);          

    %**********************************************************************
    %***
    %***     Generate graph
    %***
    %**********************************************************************
    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    figure(1);
    clf;
    
    xmin = 1860;
    xmax = 1970;   
    %--------------------------------------------------------%
    % Panel (a)
    subplot(4,4,1)
    plot(1909:1:1970,rgnp,'-k','LineWidth',0.75);
    title('RGNP','FontSize',8);
    axis manual
    xlim([xmin xmax]) 
    set(gca,'xtick',[1880 1920 1960]);
    box off
    set(gca,'FontSize',8);
        
    %--------------------------------------------------------%
    % Panel (b)
    subplot(4,4,2)
    plot(1909:1:xmax,gnp,'-k','LineWidth',0.75);
    title('GNP','FontSize',8);
    axis manual
    xlim([xmin xmax]) 
    set(gca,'xtick',[1880 1920 1960]);
    box off
    set(gca,'FontSize',8);

    %--------------------------------------------------------%
    % Panel (c)
    subplot(4,4,3)
    plot(1909:1:xmax,pcrgnp,'-k','LineWidth',0.75);
    title('PCRGNP','FontSize',8);
    axis manual
    xlim([xmin xmax]) 
    set(gca,'xtick',[1880 1920 1960]);
    box off
    set(gca,'FontSize',8); 

    %--------------------------------------------------------%
    % Panel (d)
    subplot(4,4,4)
    plot(xmin:1:xmax,ip,'-k','LineWidth',0.75);
    title('IP','FontSize',8);
    axis manual
    xlim([xmin xmax]) 
    set(gca,'xtick',[1880 1920 1960]);
    box off
    set(gca,'FontSize',8);

    %--------------------------------------------------------%
    % Panel (e)
    subplot(4,4,5)
    plot(1890:1:xmax,emp,'-k','LineWidth',0.75);
    title('Employment','FontSize',8);
    axis manual
    xlim([xmin xmax]) 
    set(gca,'xtick',[1880 1920 1960]);
    box off
    set(gca,'FontSize',8);

    %--------------------------------------------------------%
    % Panel (f)
    subplot(4,4,6)
    plot(1890:1:xmax,un,'-k','LineWidth',0.75);
    title('Unemployment','FontSize',8);
    axis manual
    xlim([xmin xmax]) 
    set(gca,'xtick',[1880 1920 1960]);
    box off
    set(gca,'FontSize',8);

    %--------------------------------------------------------%
    % Panel (g)
    subplot(4,4,7)
    plot(1889:1:xmax,prgnp,'-k','LineWidth',0.75);
    title('PRGNP','FontSize',8);
    axis manual
    xlim([xmin xmax]) 
    set(gca,'xtick',[1880 1920 1960]);
    box off
    set(gca,'FontSize',8);

    %--------------------------------------------------------%
    % Panel (h)
    subplot(4,4,8)
    plot(1860:1:xmax,cpi,'-k','LineWidth',0.75);
    title('CPI','FontSize',8);
    axis manual
    xlim([xmin xmax]) 
    set(gca,'xtick',[1880 1920 1960]);
    box off
    set(gca,'FontSize',8);

    %--------------------------------------------------------%
    % Panel (i)
    subplot(4,4,9)
    plot(1900:1:xmax,wg,'-k','LineWidth',0.75);
    title('Wages','FontSize',8);
    axis manual
    xlim([xmin xmax]) 
    set(gca,'xtick',[1870 1920 1960]);
    box off
    set(gca,'FontSize',8);
    
    %--------------------------------------------------------%
    % Panel (j)
    subplot(4,4,10)
    plot(1900:1:xmax,rwg,'-k','LineWidth',0.75);
    title('Real Wages','FontSize',8);
    axis manual
    xlim([xmin xmax]) 
    set(gca,'xtick',[1880 1920 1960]);
    box off
    set(gca,'FontSize',8);
    
    %--------------------------------------------------------%
    % Panel (k)
    subplot(4,4,11)
    plot(1889:1:xmax,m,'-k','LineWidth',0.75);
    title('Money','FontSize',8);
    axis manual
    xlim([xmin xmax]) 
    set(gca,'xtick',[1880 1920 1960]);
    box off
    set(gca,'FontSize',8);
    
    %--------------------------------------------------------%
    % Panel (l)
    subplot(4,4,12)
    plot(1869:1:xmax,vel,'-k','LineWidth',0.75);
    title('Velocity','FontSize',8);
    axis manual
    xlim([xmin xmax]) 
    set(gca,'xtick',[1880 1920 1960]);
    box off
    set(gca,'FontSize',8);
    
    %--------------------------------------------------------%
    % Panel (m)
    subplot(4,4,13)
    plot(1900:1:xmax,bnd,'-k','LineWidth',0.75);
    title('Bond Yield','FontSize',8);
    axis manual
    xlim([xmin xmax]) 
    set(gca,'xtick',[1880 1920 1960]);
    box off
    set(gca,'FontSize',8);
     
    %--------------------------------------------------------%
    % Panel (n)
    subplot(4,4,14)
    plot(1871:1:xmax,sp500,'-k','LineWidth',0.75);
    title('$S\&P500$','FontSize',8);
    axis manual
    xlim([xmin xmax]) 
    set(gca,'xtick',[1880 1920 1960]);
    box off
    set(gca,'FontSize',8);

    %--------------------------------------------------------%
    % Panel (o)
    subplot(4,4,15)
    plot(xmin:1:xmax,y_dtrend,'-k','LineWidth',0.75);
    title('Deterministic Trend','FontSize',8);
    axis manual
    xlim([xmin xmax]) 
    set(gca,'xtick',[1880 1920 1960]);
    box off
    set(gca,'FontSize',8);
    
    %--------------------------------------------------------%
    % Panel (p)
    subplot(4,4,16)
    plot(xmin:1:xmax,y_strend,'-k','LineWidth',0.75);
    title('Stochastic Trend','FontSize',8);
    axis manual
    xlim([xmin xmax]) 
    set(gca,'xtick',[1880 1920 1960]);
    box off
    set(gca,'FontSize',8);
    
    % Print the tex file to the relevant directory
    %laprint(1,'nelsonplosser','options','factory','keepfontprops', 'on');

    % Compute autocorrelations for the first 6 lags
    lags = 6;
    % Nelson and Plosser (1982) Table 2, p.147
    format short
    disp('Autocorrelation functions')
    disp('Nominal GNP')
    disp('     1        2         3         4         5          6' )
    disp( acf(gnp,lags)' );
    disp('')
    disp('Per capita GNP')
    disp('     1        2         3         4         5          6' )
    disp( acf(pcrgnp,lags)' );
    disp('')
    disp('Industrial production')
    disp('     1        2         3         4         5          6' )
    disp( acf(ip,lags)' );
    disp('')
    disp('Unemployment')
    disp('     1        2         3         4         5          6' )
    disp( acf(emp,lags)' );
    disp('')
    disp('Unemployment rate')
    disp('     1        2         3         4         5          6' )
    disp( acf(un,lags)' );
    disp('')
    disp('GNP deflator')
    disp('     1        2         3         4         5          6' )
    disp( acf(prgnp,lags)' );
    disp('')
    disp('Consumer price index')
    disp('     1        2         3         4         5          6' )
    disp( acf(cpi,lags)' );
    disp('')
    disp('Wages')
    disp('     1        2         3         4         5          6' )
    disp( acf(wg,lags)' );
    disp('')
    disp('Real wages')
    disp('     1        2         3         4         5          6' )
    disp( acf(rwg,lags)' );
    disp('')
    disp('Money stock')
    disp('     1        2         3         4         5          6' )
    disp( acf(m,lags)' );
    disp('')
    disp('Velocity')
    disp('     1        2         3         4         5          6' )
    disp( acf(vel,lags)' );
    disp('')
    disp('Bond yield')
    disp('     1        2         3         4         5          6' )
    disp( acf(bnd,lags)' );
    disp('')
    disp('Stock price')
    disp('     1        2         3         4         5          6' )
    disp( acf(sp500,lags)' );
    disp('')
    disp('Deterministic trend')
    disp('     1        2         3         4         5          6' )
    disp( acf(y_dtrend,lags)' );
    disp('')
    disp('Stochastic trend')
    disp('     1        2         3         4         5          6' )
    disp( acf(y_strend,lags)' );
    disp('')

    % Estimate the stochastic trend model 
    % an AR(1) model with a constant and no time tend 
    format short
    [b,se] = ar(rgnp,0);
    disp('Real GNP')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['AR1 parameter    ',num2str(b(2)),'  (',num2str(se(2)),')'])
   
    [b,se] = ar(gnp,0);
    disp(' ')
    disp('Nominal GNP')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['AR1 parameter    ',num2str(b(2)),'  (',num2str(se(2)),')'])
    [b,se] = ar(pcrgnp,0);
    disp(' ')
    disp('Per capita GNP')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['AR1 parameter    ',num2str(b(2)),'  (',num2str(se(2)),')'])
    [b,se] = ar(ip,0);
    disp(' ')
    disp('Industrial production')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['AR1 parameter    ',num2str(b(2)),'  (',num2str(se(2)),')'])
    [b,se] = ar(emp,0);
    disp(' ')
    disp('Employment')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['AR1 parameter    ',num2str(b(2)),'  (',num2str(se(2)),')'])
    [b,se] = ar(un,0);
    disp(' ')
    disp('Unemployment rate')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['AR1 parameter    ',num2str(b(2)),'  (',num2str(se(2)),')'])
    [b,se] = ar(prgnp,0);
    disp(' ')
    disp('GNP deflator')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['AR1 parameter    ',num2str(b(2)),'  (',num2str(se(2)),')'])
    [b,se] = ar(cpi,0);
    disp(' ')
    disp('Consumer price index')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['AR1 parameter    ',num2str(b(2)),'  (',num2str(se(2)),')'])
    [b,se] = ar(wg,0);
    disp(' ')
    disp('Wages')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['AR1 parameter    ',num2str(b(2)),'  (',num2str(se(2)),')'])
    [b,se] = ar(rwg,0);
    disp(' ')
    disp('Real wages')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['AR1 parameter    ',num2str(b(2)),'  (',num2str(se(2)),')'])
    [b,se] = ar(m,0);
    disp(' ')
    disp('Money stock')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['AR1 parameter    ',num2str(b(2)),'  (',num2str(se(2)),')'])
    [b,se] = ar(vel,0);
    disp(' ')
    disp('Velocity')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['AR1 parameter    ',num2str(b(2)),'  (',num2str(se(2)),')'])
    [b,se] = ar(bnd,0);
    disp(' ')
    disp('Bond yield')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['AR1 parameter    ',num2str(b(2)),'  (',num2str(se(2)),')'])
    [b,se] = ar(sp500,0);
    disp(' ')
    disp('Stock price')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['AR1 parameter    ',num2str(b(2)),'  (',num2str(se(2)),')'])
    [b,se] = ar(y_dtrend,0);
    disp(' ')
    disp('Deterministic trend')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['AR1 parameter    ',num2str(b(2)),'  (',num2str(se(2)),')'])
    [b,se] = ar(y_strend,0);
    disp(' ')
    disp('Stochastic trend')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['AR1 parameter    ',num2str(b(2)),'  (',num2str(se(2)),')'])

    % Estimate the combined model: 
    % an AR(1) model with a constant and a time tend 
    [b,se] = ar(rgnp,1);
    disp('Real GNP')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    disp(['AR1 parameter    ',num2str(b(3)),'  (',num2str(se(3)),')']) 
    [b,se] = ar(gnp,1);
    disp(' ')
    disp('Nominal GNP')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    disp(['AR1 parameter    ',num2str(b(3)),'  (',num2str(se(3)),')'])
    [b,se] = ar(pcrgnp,1);
    disp(' ')
    disp('Per capita GNP')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    disp(['AR1 parameter    ',num2str(b(3)),'  (',num2str(se(3)),')'])
    [b,se] = ar(ip,1);
    disp(' ')
    disp('Industrial production')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    disp(['AR1 parameter    ',num2str(b(3)),'  (',num2str(se(3)),')'])
    [b,se] = ar(emp,1);
    disp(' ')
    disp('Employment')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    disp(['AR1 parameter    ',num2str(b(3)),'  (',num2str(se(3)),')'])
    [b,se] = ar(un,1);
    disp(' ')
    disp('Unemployment rate')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    disp(['AR1 parameter    ',num2str(b(3)),'  (',num2str(se(3)),')'])
    [b,se] = ar(prgnp,1);
    disp(' ')
    disp('GNP deflator')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    disp(['AR1 parameter    ',num2str(b(3)),'  (',num2str(se(3)),')'])
    [b,se] = ar(cpi,1);
    disp(' ')
    disp('Consumer price index')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    disp(['AR1 parameter    ',num2str(b(3)),'  (',num2str(se(3)),')'])
    [b,se] = ar(wg,1);
    disp(' ')
    disp('Wages')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    disp(['AR1 parameter    ',num2str(b(3)),'  (',num2str(se(3)),')'])
    [b,se] = ar(rwg,1);
    disp(' ')
    disp('Real wages')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    disp(['AR1 parameter    ',num2str(b(3)),'  (',num2str(se(3)),')'])
    [b,se] = ar(m,1);
    disp(' ')
    disp('Money stock')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    disp(['AR1 parameter    ',num2str(b(3)),'  (',num2str(se(3)),')'])
    [b,se] = ar(vel,1);
    disp(' ')
    disp('Velocity')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    disp(['AR1 parameter    ',num2str(b(3)),'  (',num2str(se(3)),')'])
    [b,se] = ar(bnd,1);
    disp(' ')
    disp('Bond yield')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    disp(['AR1 parameter    ',num2str(b(3)),'  (',num2str(se(3)),')'])
    [b,se] = ar(sp500,1);
    disp(' ')
    disp('Stock price')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    disp(['AR1 parameter    ',num2str(b(3)),'  (',num2str(se(3)),')'])
    [b,se] = ar(y_dtrend,1);
    disp(' ')
    disp('Deterministic trend')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    disp(['AR1 parameter    ',num2str(b(3)),'  (',num2str(se(3)),')'])
    [b,se] = ar(y_strend,1);
    disp(' ')
    disp('Stochastic trend')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    disp(['AR1 parameter    ',num2str(b(3)),'  (',num2str(se(3)),')'])
 
    % Estimate the trend model: 
    % a model with a constant and a time tend 
   [b,se] = ar(rgnp,2);
    disp('Real GNP')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])       
    [b,se] = ar(gnp,2);
    disp(' ')
    disp('Nominal GNP')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    [b,se] = ar(pcrgnp,2);
    disp(' ')
    disp('Per capita GNP')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    [b,se] = ar(ip,2);
    disp(' ')
    disp('Industrial production')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    [b,se] = ar(emp,2);
    disp(' ')
    disp('Employment')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    [b,se] = ar(un,2);
    disp(' ')
    disp('Unemployment rate')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    [b,se] = ar(prgnp,2);
    disp(' ')
    disp('GNP deflator')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    [b,se] = ar(cpi,2);
    disp(' ')
    disp('Consumer price index')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    [b,se] = ar(wg,2);
    disp(' ')
    disp('Wages')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    [b,se] = ar(rwg,2);
    disp(' ')
    disp('Real wages')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    [b,se] = ar(m,2);
    disp(' ')
    disp('Money stock')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    [b,se] = ar(vel,2);
    disp(' ')
    disp('Velocity')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    [b,se] = ar(bnd,2);
    disp(' ')
    disp('Bond yield')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    [b,se] = ar(sp500,2);
    disp(' ')
    disp('Stock price')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    [b,se] = ar(y_dtrend,2);
    disp(' ')
    disp('Deterministic trend')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    
    [b,se] = ar(y_strend,2);
    disp(' ')
    disp('Stochastic trend')
    disp(['constant         ',num2str(b(1)),'  (',num2str(se(1)),')'])
    disp(['trend            ',num2str(b(2)),'  (',num2str(se(2)),')'])    

end
%--------------------------- Functions  ----------------------------------
% 
%--------------------------------------------------------------------------
% ACF function based on running a sequence of OLS regressions  
%--------------------------------------------------------------------------
function acorr = acf(y,lags)

    acorr = zeros(lags,1);
    t     = length(y);

    for i = 1:lags
        
        y0 = trimr(y,i,0);
        x  = [ ones(t-i,1) trimr(y,0,i)];
        b  = inv(x'*x)*x'*y0;
        
        acorr(i) = b(2);
    end
end
%--------------------------------------------------------------------------
% Regression estimates of trend models
%   ntrend = 0 (constant and lag)
%   ntrend = 1 (constant, time trend and lag)
%   ntrend = 2 (constant and time trend)
%--------------------------------------------------------------------------
function [ b,se ] = ar(y,ntrend)

    t = length(y);
    
    if ntrend == 0.0
        x = [ ones(t-1,1)   trimr(y,0,1) ];                                
        y = trimr(y,1,0);                                                    
    elseif ntrend == 1.0
        x = [ ones(t-1,1)   seqa(0,1,t-1)'/100   trimr(y,0,1) ];      
        y = trimr(y,1,0);                                                      
    else
        x = [ ones(t,1)   seqa(0,1,t)'/100 ];                          
        y = trimr(y,0,0);                                                           
    end

    b = x\y;                                    
    e = y - x*b;                                 
    s2 = e'*e/length(e);                            
    se = sqrt( diag( s2*inv(x'*x) ) );          

end
