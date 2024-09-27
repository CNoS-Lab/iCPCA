function [fig] = plotScree(CPCA_struct)

fig = figure('Position',[0 0 1280 720],'visible','on'); % Scree plot
D = CPCA_struct.DGH;
D = diag(D).^2;
D = D./sum(D)*100;
k = length(D);
plot(1:k,D(1:k),'-O');
xlabel('No. of Components');
ylabel('Percentage of variance');
title(sprintf('Scree Plot'));
set(findall(gcf,'-property','FontSize'),'FontSize',15);

end