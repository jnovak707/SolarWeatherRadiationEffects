function [pd1,pd2,pd3] = createFit(NOAA_PFU_Array,test)
%CREATEFIT    Create plot of datasets and fits
%   [PD1,PD2,PD3] = CREATEFIT(NOAA_PFU_ARRAY,TEST)
%   Creates a plot, similar to the plot in the main distribution fitter
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with distributionFitter
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  3
%   Number of fits:  3
%
%   See also FITDIST.

% This function was automatically generated on 10-Apr-2022 19:31:59

% Output fitted probablility distributions: PD1,PD2,PD3

% Data from dataset "NOAA_PFU_Array data":
%    Y = NOAA_PFU_Array

% Data from dataset "test data":
%    Y = test

% Data from dataset "All SPE PFU (log10 Scale)":
%    Y = test

% Force all inputs to be column vectors
NOAA_PFU_Array = NOAA_PFU_Array(:);
test = test(:);

% Prepare figure
clf;
hold on;
LegHandles = []; LegText = {};
set(gcf,'DefaultAxesFontName','Times New Roman','DefaultAxesFontWeight','bold');
set(gcf,'DefaultTextFontName','Times New Roman','DefaultAxesFontWeight','bold');
set(gcf,'DefaultTextColor','black')

% --- Plot data originally in dataset "NOAA_PFU_Array data"
% This dataset does not appear on the plot

% --- Plot data originally in dataset "test data"
% This dataset does not appear on the plot

% --- Plot data originally in dataset "All SPE PFU (log10 Scale)"
[CdfF,CdfX] = ecdf(test,'Function','cdf');  % compute empirical cdf
BinInfo.rule = 5;
BinInfo.width = 1;
BinInfo.placementRule = 1;
[~,BinEdge] = internal.stats.histbins(test,[],[],BinInfo,CdfF,CdfX);
[BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
hLine = bar(BinCenter,BinHeight,'hist');
set(hLine,'FaceColor','none','EdgeColor',[0 0 0],...
    'LineStyle','-', 'LineWidth',1);
xlabel('Solar Radiation Storm Level (S1-S5)','FontName','Times','FontWeight','bold', 'FontSize',26');
ylabel('Density','FontName','Times','FontWeight','bold', 'FontSize',26');
LegHandles(end+1) = hLine;
LegText{end+1} = 'All SPE PFU (log10 Scale)';

% Create grid where function will be computed
XLim = get(gca,'XLim');
XLim = XLim + [-1 1] * 0.01 * diff(XLim);
XGrid = linspace(XLim(1),XLim(2),100);
a= get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',26)



% --- Create fit "Non_Parametric"
pd1 = fitdist(NOAA_PFU_Array,'kernel','kernel','normal','support',[9, 110000],'width',100);
% This fit does not appear on the plot

% --- Create fit "fit 2"
pd2 = fitdist(test,'kernel','kernel','normal','support','unbounded');
% This fit does not appear on the plot

% --- Create fit "Non-Parametric KDE (Normal)"
pd3 = fitdist(test,'kernel','kernel','normal','support','unbounded');
YPlot = cdf(pd3,XGrid);
hLine = plot(XGrid,YPlot,'Color',[0.9 0 0],...
    'LineStyle','-', 'LineWidth',2,...
    'Marker','none', 'MarkerSize',6);
LegHandles(end+1) = hLine;
LegText{end+1} = 'Non-Parametric KDE (Normal)';

% Adjust figureXGrid
box on;
grid on
hold off;

% Create legend from accumulated handles and labels
hLegend = legend(LegHandles,LegText,'Orientation', 'vertical', 'FontSize', 26,'FontWeight','bold', 'Location', 'northeast');
set(hLegend,'Interpreter','none');


xlim([1 5.1]);

