function [ax2] = add_xaxis(ax1,activated,xlimval,xtickval,xticklab)
%ADD_XAXIS Add a second x-axis on top of the plot with the specified
%tickvalues.
%
%   INPUTS:
%   ax1: handle of the existing axis.
%   activated (string): 'on' (the new axis is activated) or 'off' (the original
%   axis is activated).
%   xlimval (numerical array, optional): limits of the axis.
%   xtickval (numerical array, optional): values of the ticks.
%   xticklab (cell array, optional): labels of the thicks.
%
%   OUTPUTS:
%   ax2: handle of the new axis
%
% T. Doda, 06.01.2022
%%

axpos=get(ax1,'Position');
ax2=axes('Position',axpos); hold on;
ax2.YLimMode='manual';
if strcmp(ax1.YAxisLocation,'left')
    ax2.YAxisLocation='right';
else
    ax2.YAxisLocation='left';
end

if nargin>=4
set(ax2,'xlim',xlimval,'xtick',xtickval,'ylim',ax1.YLim,'ytick',ax1.YTick,...
    'yticklabel','','XAxisLocation', 'top',...
    'box','off','fontsize',ax1.FontSize)
else
    set(ax2,'ylim',ax1.YLim,'ytick',ax1.YTick,'yticklabel','',...
    'XAxisLocation', 'top',...
    'box','off','fontsize',ax1.FontSize)
end
set(ax1,'box','off');
if nargin==5
    ax2.XTickLabel=xticklab;
end

if strcmp(activated,'off')
    set(ax1,'color','none')
    axes(ax1)
else
    set(ax2,'color','none');
end
end

