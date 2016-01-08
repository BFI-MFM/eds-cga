function MakeAllFigures(NumberOfPlots,Caption,FigureProperties,Info)

% Copyright (C) 2005-2009 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

global M_ options_

FigHandle = figure('Name',FigureProperties.Name);  

NAMES = cell(NumberOfPlots,1);
if options_.TeX
    TeXNAMES = cell(NumberOfPlots,1); 
end

if NumberOfPlots == 9
    nr = 3;
    nc = 3;
elseif NumberOfPlots == 8
    nr = 3;
    nc = 3;
elseif NumberOfPlots == 7
    nr = 3;
    nc = 3;
elseif NumberOfPlots == 6
    nr = 2;
    nc = 3;
elseif NumberOfPlots == 5
    nr = 3;
    nc = 2;
elseif NumberOfPlots == 4
    nr = 2;
    nc = 2;
elseif NumberOfPlots == 3
    nr = 2;
    nc = 2;
elseif NumberOfPlots == 2
    nr = 1;
    nc = 2;
elseif NumberOfPlots == 1
    nr = 1;
    nc = 1;
end  

for plt = 1:NumberOfPlots
    eval(['NumberOfCurves = Info.Box' int2str(plt) '.Number;'])
    NumberOfObservations = zeros(2,1);
    x = cell(NumberOfCurves,1);
    y = cell(NumberOfCurves,1);
    PltType = cell(NumberofCurves,1);
    top = NaN(NumberOfCurves,1);
    bottom = NaN(NumberOfCurves,1);
    binf = NaN(NumberOfCurves,1);
    bsup = NaN(NumberOfCurves,1);
    for curve = 1:NumberOfCurves
        eval(['x{' curve '} = Info.Box' int2str(plt) '.Curve' int2str(curve) '.xdata;'])
        eval(['y{' curve '} = Info.Box' int2str(plt) '.Curve' int2str(curve) '.ydata;'])
        eval(['name = Info.Box' int2str(plt) '.Curve' int2str(curve) '.variablename;'])
        eval(['PltType{' curve  '} = Info.Box' int2str(plt) '.Curve' int2str(curve) '.type']);
        if length(x{curve})-length(y{curve})
            disp('MakeFigure :: The number of observations in x doesn''t match with ')
            disp(['the number of observation in y for ' name ])
            return
        end
        if Info.PlotProperties.CutTop
            top(curve) = max(y{curve});
        else Info.PlotProperties.CutBottom
            bottom(curve) = min(y{curve});
        end
        binf(curve) = min(x{curve});
        bsup(curve) = max(x{curve});
    end
    ymax = max(top);
    ymin = min(bottom);
    xmin = min(binf);
    xmax = max(bsup);
    if isnan(ymin(plt))
        ymin = 0;
    end
    eval(['NAMES{' int2str(plt) '} = Info.Box' int2str(plt) '.name;'])
    if options_.TeX
        eval(['TeXNAMES{' int2str(plt) '} = Info.Box' int2str(plt) '.texname;'])
    end
    subplot(nr,nc,plt)
    hold on
    for curve = 1:NumberOfCurves
        hh = plot(x{curve},y{curve});
        if strcmpi(PltType{curve},'PriorDensity')
            set(hh,'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',2)
            %
            %
        elseif strcmpi(PltType{curve},'DensityEstimate')
            set(hh,'Color','k','LineStyle','-','LineWidth',2)
            %
            %
        elseif strcmpi(PltType{curve},'ModeEstimate')
            set(hh,'Color','g','LineStyle','--','LineWidth',2)
            %
            %
        elseif strcmpi(PltType{curve},'SmoothVariable')
            set(hh,'Color','k','LineStyle','-','LineWidth',2)
            %
            %
        elseif strcmpi(PltType{curve},'Deciles')
            set(hh,'Color','g','LineStyle','-','LineWidth',1)
            %
            %
        elseif strcmpi(PltType{curve},'Forecasts')
            set(hh,'Color','','LineStyle','-','LineWidth',2)
            %
            %
        elseif strcmpi(PltType{curve},'ForecastsHPD')
            set(hh,'Color','k','LineStyle','-','LineWidth',1)
            %
            %
        elseif strcmpi(PltType{curve},'ForecastsDeciles')
            set(hh,'Color','g','LineStyle','-','LineWidth',1)
            %
            %
        elseif strcmpi(PltType{curve},'DiagnosticWithin')
            set(hh,'Color','b','LineStyle','-','LineWidth',2)
            %
            %
        elseif strcmpi(PltType{curve},'DiagnosticPooled')
            set(hh,'Color','r','LineStyle','-','LineWidth',2)
            %
            %
        end  
    end
    axis([xmin xmax ymin ymax])
    title(NAMES{plt})
    drawnow
    hold off
end

if Info.SaveFormat.Eps
    if isempty(Info.SaveFormat.Name)
        eval(['print -depsc2 ' M_.fname Info.SaveFormat.GenericName int2str(Info.SaveFormat.Number) '.eps']);
    else
        eval(['print -depsc2 ' M_.fname Info.SaveFormat.GenericName Info.SaveFormat.Name '.eps']);  
    end
end
if Info.SaveFormat.Pdf && ~exist('OCTAVE_VERSION')
    if isempty(Info.SaveFormat.Name)
        eval(['print -dpdf ' M_.fname Info.SaveFormat.GenericName int2str(Info.SaveFormat.Number)]);
    else
        eval(['print -dpdf ' M_.fname Info.SaveFormat.GenericName Info.SaveFormat.Name]);  
    end
end
if Info.SaveFormat.Fig && ~exist('OCTAVE_VERSION')
    if isempty(Info.SaveFormat.Name)
        saveas(FigHandle,[M_.fname Info.SaveFormat.GenericName int2str(Info.SaveFormat.Number) '.fig']);
    else
        saveas(FigHandle,[M_.fname Info.SaveFormat.GenericName Info.SaveFormat.Name '.fig']);
    end
end