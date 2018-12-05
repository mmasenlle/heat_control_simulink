%% Plot salidas y acciones
close all
figure
h(1) = subplot(4,1,1); plot(t_sim,x_sim(:,params.op)),legend('T1'),grid
h(2) = subplot(4,1,2); plot(t_sim,q_sim),legend('Pot'),grid
h(3) = subplot(4,1,3); plot(t_sim,x_sim(:,params.op2)),legend('T2'),grid
h(4) = subplot(4,1,4); plot(t_sim,qv_sim),legend('vel'),grid
linkaxes(h,'x')