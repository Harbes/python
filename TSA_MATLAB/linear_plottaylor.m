%=========================================================================
%
%   Program to plot data from Taylor rules
%
%=========================================================================

function linear_plottaylor()


    clear;
    clc;
    
    %Read data from timeseries obj
	load taylor;
        
    % Storing data on arrays.
	infl  = taylor.Data(:,1);
	ygap  = taylor.Data(:,2);
	ffr   = taylor.Data(:,3);
    
    
    % Get the time variable for plotting purposes
    dates = taylor.Time;

    %**********************************************************************
    %***
    %***     Generate graph
    %***
    %**********************************************************************

    % Switch off TeX interpreter and clear figure
    set(0,'defaulttextinterpreter','none');
    f1=figure(1);
    clf;

    plot(dates,infl,'-k',dates,ffr,'--k',dates,ygap,':k','LineWidth',0.75);
    ylabel('Percent')
    datetick('x','yyyy');
    %legend('inflation','fed funds rate','ouput gap','Location','NorthWest')
    %legend2latex( f1 )
    axis tight
    box off


    % Print the tex file to the relevant directory
    laprint(1,'taylor','options','factory');


end