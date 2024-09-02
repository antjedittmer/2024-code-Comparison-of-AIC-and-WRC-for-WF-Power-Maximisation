function axislesstight(ax)
if ~nargin
    ax = gca;
end
ylim = ax.YLim;
dYlim = ylim(2) -ylim(1);
ax.YLim = [ylim(1) - 0.1*dYlim, ylim(2) + 0.1*dYlim];