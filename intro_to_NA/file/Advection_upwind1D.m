%% Input Parameter
c       = 1.33;
x_left  = -3;
x_right = 5*10;
x_num   = 5*50;
x       = linspace(x_left, x_right, x_num);
dx      = (x_right-x_left)/(x_num-1);
% phi_init= 0*x; phi(10:15) = 1.0; % Potition data
phi_init= exp(-x.^2);     

dt      = 0.2;
t_max   = 15;
t       =  0;

%% Plot
fig = figure;

% %Normal plot
% p = plot(x, phi_init); % Initial plot
% xlim([x_left,x_right]);ylim([-0.1,1.1]);
% name = title(['t = ',num2str(t),'[s]']);

% % compare plot
% ax for exact solution
ax_exact = subplot(2,1,1);
name_exact = title(ax_exact,['解析解:t = ',num2str(t,'%04.1f'),'[s]']);
ax_exact.XLim = [x_left, x_right];
ax_exact.YLim = [-0.1, 1.1];

% ax for simulation result
ax_sim   = subplot(2,1,2);
name_sim = title(ax_sim,['数値解:t = ',num2str(t,'%04.1f'),'[s]']);
ax_sim.XLim = [x_left, x_right];
ax_sim.YLim = [-0.1, 1.1];

p_exact = plot(ax_exact, x, exp(-(x-c*t).^2));
p_sim   = plot(ax_sim, x, phi_init);

%% Simulation
phi_n1 = phi_init;
phi_n  = phi_init;
h = c*dt/dx;
% stop
for t = 0:dt:t_max
    for i = 2:x_num-1
        phi_n1(i) = (1-h)*phi_n(i)+h*phi_n(i-1);
    end
    % update boundary
    phi_n1(1)   = phi_n1(2);
    phi_n1(end) = phi_n1(end-1);
    
    % update plot
%     p.YData  = phi_n1;
%     name = title(['t = ',num2str(t),'[s]']);
    name_exact = title(ax_exact,['解析解:t = ',num2str(t,'%04.1f'),'[s]']);
    name_sim = title(ax_sim,['数値解:t = ',num2str(t,'%04.1f'),'[s]']);
    p_sim.YData = phi_n1;
    p_exact.YData = exp(-(x-c*t).^2);
    ax_sim.YLim = [-0.1, 1.1];
    ax_sim.XLim = [x_left, x_right];
    ax_exact.YLim = [-0.1, 1.1];
    ax_exact.XLim = [x_left, x_right];
    
    drawnow;
    % update n
    
    phi_n = phi_n1;
    % save figure
%     filename = sprintf('result/advection/upwind/plotUPWIND_%04.1f.png', t); % file name 
    filename = sprintf('result/advection/compare/plotUPWIND_compareV%04.1f.png', t); % file name 
%     filename = sprintf('result/advection/upwind/plot_vibration%04.1f.png', t); % file name 
%     saveas(fig, filename); % save figure as png    

end