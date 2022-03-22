%% Input Parameters
c = 1;
t = 0;
x = linspace(-3,10, 200);
fig = figure;
p = plot(x, exp(-(x-c*t).^2)); % Initial plot
xlim([-3,10]);ylim([-0.1,1.1]);
name = title(['t = ',num2str(t),'[s]']);


for t = 0:0.2:15 % time data
    name = title(['t = ',num2str(t,'%04.1f'),'[s]']);
    p.YData = exp(-(x-c*t).^2);
    drawnow;
    filename = sprintf('plotEXACT%04.1f.png', t); % file name 
    %saveas(fig, filename); % save figure as png
end