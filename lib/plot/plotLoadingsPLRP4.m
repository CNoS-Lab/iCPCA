    function [ax] = plotLoadingsPLRP4(X,PLRP,xTicks,yTicks,pValues,p_thresh,cutoff)

h = imagesc(X);
xticklabels(xTicks);
yticklabels(yTicks);

ax = gca;
ax.XTick = 1:size(X,2);
ax.XTickLabel = xTicks;
ax.YTick = 1:size(X,1);
ax.YTickLabel = yTicks;     

[rows,cols] = size(X);

for i = 1:rows
    for j = 1:cols
        if pValues(i,j) <= p_thresh(3)
            pStr = '***';
        elseif pValues(i,j) <= p_thresh(2)
            pStr = '**';
        elseif pValues(i,j) <= p_thresh(1)
            pStr = '*';
        else
            pStr = '';
        end

        if PLRP(i,j) >= cutoff
            xStr = '*';
        else
            xStr = '';
        end

        textHandles(j,i) = text(j,i,sprintf('%.2f\n(%.2f%s, \n%.3f%s)',X(i,j),PLRP(i,j),xStr,pValues(i,j),pStr),...
                'horizontalAlignment','center');
        

    end
end

ax = gca;

end