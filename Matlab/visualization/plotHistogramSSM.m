%================================================
%
%  function plotHistogramSSM
%
%=================================================
%  INPUT:
%  data     = vector of data for which to plot histogram
%  nBind    = number of bins to use (fun optBINS.m finds optimal number of bins)
%  bBarPlot = 1 plot histogram as bar, 0 as points
%  col      = color of the histgram
%  lw       = line widht
%  ms       = marker size
%  mar      = type of the marker for the point plot, pass as a string


function plotHistogramSSM(varargin)

data        = varargin{1};
nBins       = varargin{2};
bBarPlot    = varargin{3};
col         = varargin{4};
opacity     = varargin{5};
lw          = varargin{6};

if (nargin > 6)
    mar = varargin{7};
    ms  = varargin{8};
end;

if(bBarPlot)
    histogram(data,nBins,'Normalization','pdf','EdgeColor',col,'FaceColor',col,'FaceAlpha',opacity,'LineWidth',lw);
else
    [values,edges] = histcounts(data,nBins,'Normalization','pdf');
    halfWidth = 0.5*(edges(2) - edges(1));
    binCenters     = linspace( edges(1) + halfWidth, edges(end)-halfWidth, nBins);
    plot(binCenters,values, mar, 'MarkerSize',ms,'LineWidth',lw,'Color',col);
end;



