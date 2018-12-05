tend=10;

% linealizaciones
model = 'lin_2x2';
opspec = operspec(model)
opspec.Inputs(1).u = params.Pot_nom;
opspec.Inputs(2).u = params.v_nom;
opspec.Inputs(1).Known = true;
opspec.Inputs(2).Known = true;
[op,opreport] = findop(model,opspec);
linsys = linearize(model,op);

if 1
w = logspace(-2,2,100);
Pf = freqresp(linsys,w);

rgas=[];
for ii=1:100
P_frd = frd(Pf,w);
rga=Pf(:,:,ii).*inv(Pf(:,:,ii)).';
rga=real(rga);

rgas=[rgas; rga(1,1), rga(1,2)];
end

return

% plot([1,2,3,4],[rga(1,1),rga(1,2),rga(2,1),rga(2,2)])
% fname=sprintf('img/rga_xi_%d__xj_%d.png', xi_mm, xj_mm);
% saveas(gcf,fname);
% close all
end

%%
% equilibrio, mirar a ojo enfriamiento, dependera de conveccion (params.h)
if 0
    figure,hold on
    m=min(min(op.state.x)); M=max(max(op.state.x))*1.2;
    plot(nodes,op.state.x)
    plot([params.input_x params.input_x],[m M],'m--');
    plot(nodes([params.op params.op]),[m M],'y--');
    plot(nodes([params.op2 params.op2]),[m M],'c--');
    hold off
end

% Validaci�n del modelo lineal MIMO
if 1
tend = 50; 
%[ylin,tlin] = step(linsys,tend);

%sim_init
t_step = 0;
delta_pot = 500; delta_vel = 0;
sim open_loop;

t_sim1=[t_sim1 t_sim];
x_sim1=[x_sim1 x_sim(:,params.op)];
t_sim3=[t_sim3 t_sim];
x_sim3=[x_sim3 x_sim(:,params.op2)];

%figure(2); clf
% h(1) = subplot(2,2,1); 
% plot(t_sim,x_sim(:,params.op),tlin,ylin(:,1,1)*delta_pot+x_sim(1,params.op))
% % ylim([2015 2035])
% h(3) = subplot(2,2,3); 
% plot(t_sim,x_sim(:,params.op2),tlin,ylin(:,2,1)*delta_pot+x_sim(1,params.op2))
% % ylim([1940 1960])
% 
% %sim_init
t_step = 0;
delta_pot = 0; delta_vel = 0.00001;
sim open_loop;

t_sim2=[t_sim2 t_sim];
x_sim2=[x_sim2 x_sim(:,params.op)];
t_sim4=[t_sim4 t_sim];
x_sim4=[x_sim4 x_sim(:,params.op2)];

% 
% h(2) = subplot(2,2,2); 
% plot(t_sim,x_sim(:,params.op),tlin,ylin(:,1,2)*delta_vel+x_sim(1,params.op))
% % ylim([2015 2035])
% h(4) = subplot(2,2,4); 
% plot(t_sim,x_sim(:,params.op2),tlin,ylin(:,2,2)*delta_vel+x_sim(1,params.op2))
% % ylim([1940 1960])
% 
% fname=sprintf('img/step_xi_%d__xj_%d.png', xi_mm, xj_mm);
% saveas(gcf,fname);
% close all
end

w = logspace(-5,1,300);
fr=freqresp(linsys(1,1),w);
fr1=[fr1 abs(squeeze(fr))];
fr=freqresp(linsys(1,2),w);
fr2=[fr2 abs(squeeze(fr))];
fr=freqresp(linsys(2,1),w);
fr3=[fr3 abs(squeeze(fr))];
fr=freqresp(linsys(2,2),w);
fr4=[fr4 abs(squeeze(fr))];
ww=[ww w'];


return

%% Respuesta frecuencial y comparativa con funciones anal�ticas

% anal�tica
v=params.v_nom;
a=params.ks/(params.rho*params.cs);
b=params.P*params.h*a/(params.ks*params.Area);
xi=nodes(params.op) - params.input_x;
w = logspace(-5,1,300);
s = j*w;

% V�lido para x>0: 
resp_pot = 1/(params.rho*params.cs)*( exp((xi/(2*a)*(v-sqrt(v^2 + 4*a*b + 4*a*s)))))./sqrt(v^2 + 4*a*b + 4*a*s);
% resp_vel = params.Pot_nom*resp_pot.*((v-sqrt(v^2 + 4*a*b + 4*a*s))./(2*a*s + 2*a*b))*.33;

resp_vel = params.Pot_nom/params.ks.*(exp((v*xi)./(2*a)).*(exp(-(xi*(v^2 + 4*a*b).^(1/2))./(2*a))...
    - exp(-(xi*(v^2 + 4*a*b + 4*a*s).^(1/2))./(2*a)) - (v*exp(-(xi*(v^2 + 4*a*b).^(1/2))./(2*a)))...
    /(v^2 + 4*a*b).^(1/2) + (v*exp(-(xi*(v^2 + 4*a*b + 4*a*s).^(1/2))./(2*a)))...
    ./(v^2 + 4*a*b + 4*a*s).^(1/2)))./(2*s);

sys_pot = frd(resp_pot, w);
sys_vel = frd(resp_vel, w);

bode_opt = bodeoptions;
bode_opt.PhaseMatching = 'on';

%figure(4)
% subplot(1,2,1); bode(linsys(1,1),sys_pot,w,bode_opt)
% subplot(1,2,2); bode(linsys(1,2),sys_vel,w,bode_opt)
% ******* Nota:
%  Cuando la conveccion es alta el bode analitico de velocidad no parte
%  con pendiente de ganancia nula por problema de calculos numericos
% legend simulink analytical

subplot(2,2,1); bodemag(linsys(1,1),w,bode_opt);
subplot(2,2,2); bodemag(linsys(1,2),w,bode_opt);
subplot(2,2,3); bodemag(linsys(2,1),w,bode_opt);
subplot(2,2,4); bodemag(linsys(2,2),w,bode_opt);

fname=sprintf('img/bode_xi_%d__xj_%d.png', xi_mm, xj_mm);
saveas(gcf,fname);
close all

