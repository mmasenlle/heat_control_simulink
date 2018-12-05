% script para diseÃ±ar control de potencia
clear all;
% para modificar parametros fisicos editar sim_init
tend=10;
sim_init;

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
w = logspace(-10,2,100);
Pf = freqresp(linsys,w);
P_frd = frd(Pf,w);
Pf(:,:,1).*inv(Pf(:,:,1)).'
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

% Validación del modelo lineal MIMO
if 1
tend = 50; 
[ylin,tlin] = step(linsys,tend);

%sim_init
t_step = 0;
delta_pot = 500; delta_vel = 0;
sim open_loop;

figure(2); clf
h(1) = subplot(2,2,1); 
plot(t_sim,x_sim(:,params.op),tlin,ylin(:,1,1)*delta_pot+x_sim(1,params.op))
% ylim([2015 2035])
h(3) = subplot(2,2,3); 
plot(t_sim,x_sim(:,params.op2),tlin,ylin(:,2,1)*delta_pot+x_sim(1,params.op2))
% ylim([1940 1960])

sim_init
t_step = 0;
delta_pot = 0; delta_vel = 0.00001;
sim open_loop;

h(2) = subplot(2,2,2); 
plot(t_sim,x_sim(:,params.op),tlin,ylin(:,1,2)*delta_vel+x_sim(1,params.op))
% ylim([2015 2035])
h(4) = subplot(2,2,4); 
plot(t_sim,x_sim(:,params.op2),tlin,ylin(:,2,2)*delta_vel+x_sim(1,params.op2))
% ylim([1940 1960])
end

%% Respuesta frecuencial y comparativa con funciones analíticas

% analítica
v=params.v_nom;
a=params.ks/(params.rho*params.cs);
b=params.P*params.h*a/(params.ks*params.Area);
xi=nodes(params.op) - params.input_x;
w = logspace(-5,1,300);
s = j*w;

% Válido para x>0: 
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

figure(4)
% subplot(1,2,1); bode(linsys(1,1),sys_pot,w,bode_opt)
% subplot(1,2,2); bode(linsys(1,2),sys_vel,w,bode_opt)
% ******* Nota:
%  Cuando la conveccion es alta el bode analitico de velocidad no parte
%  con pendiente de ganancia nula por problema de calculos numericos
% legend simulink analytical

subplot(2,2,1); bodemag(linsys(1,1),w,bode_opt); hold on
subplot(2,2,2); bodemag(linsys(1,2),w,bode_opt); hold on
subplot(2,2,3); bodemag(linsys(2,1),w,bode_opt); hold on
subplot(2,2,4); bodemag(linsys(2,2),w,bode_opt); hold on



%% Perfil de equilibrio (Rosenthal) para ver que linealizamos en una zona "manejable"

%sim_init
t_step = 0;
delta_pot = 0; delta_vel = 0;
sim open_loop;
plot_perfil

clf
plot(nodes,op.state.x)

%% Análisis RGA
% linsys=sys1d_4_1 
close all
w = logspace(-10,2,100);
Pf = freqresp(linsys,w);
P_frd = frd(Pf,w);
% figure; bode(linsys,P_frd)

for i=1%:length(w)
    RGAw(:,:,i)=Pf(:,:,i).*inv(Pf(:,:,i)).';
    RGAnum(:,:,i) = sum(sum(abs(RGAw(:,:,i) - eye(2))))
end

% figure
% subplot(2,2,1); semilogx(w,abs(squeeze(RGAw(1,1,:))))
% subplot(2,2,2); semilogx(w,abs(squeeze(RGAw(1,2,:))))
% subplot(2,2,3); semilogx(w,abs(squeeze(RGAw(2,1,:))))
% subplot(2,2,4); semilogx(w,abs(squeeze(RGAw(2,2,:))))

RGAw(:,:,1)

%% diseño controladores

close all

% Diseño para estabilidad del lazo 1,1
w=1; % Como no hay incertidumbre, basta una frecuencia (cualquiera)
W1=0.5/cos(130*pi/180/2); %MF=50º
phs=-360:1:0;
B_stab=sisobnds(1,w,W1,linsys(1,1),0,1,[],1,phs)
% Gc_pot = tf(1200,[1 0]) * tf([1/0.3 1],[1/0.2 1]);
Gc_pot = tf(1500,[1 0]);
% Gc_pot = tf(0);
w1=logspace(-3,2,1000);
lpshape(w1,B_stab,linsys(1,1),Gc_pot,phs)

% Diseño para estabilidad del lazo 2,2
B_stab=sisobnds(1,w,W1,linsys(2,2),0,1,[],1,phs)
Gc_vel = tf(1.2e-6*[1/0.3 1],[1 0]);
% Gc_vel = tf(-2e-7,[1 0])*tf([1/1.5 1],[1/2 1]);
% Gc_vel = tf(-1e-6,[1 0]);
% Gc_vel = tf(0);
lpshape(w1,B_stab,linsys(2,2),Gc_vel,phs)
w1=logspace(-3,2,1000);

%% Resultados

% Simulaciones lineales independientes
close all
tend = 200;

% Lazo 1
% Gc_pot=tf(0); % desactivar controlador de potencia
Tc1=feedback(Gc_pot*linsys(1,1), 1);
[y1,t1] = step(Tc1,tend);

% Lazo 2
% Gc_vel=tf(0); % desactivar controlador de velocidad
Tc2=feedback(Gc_vel*linsys(2,2), 1);
[y2,t2] = step(Tc2,tend);

% Simulación lineal MIMO
Gc = [Gc_pot 0;0 Gc_vel]
Tc = feedback(linsys*Gc,eye(2));
[ylin,tlin] = step(Tc,tend); %Simulación lineal

% Simulación no lineal y resultados en conjunto
Gc_pot_d=c2d(Gc_pot,params.th,'tustin');
Gc_vel_d=c2d(Gc_vel,params.th,'tustin');

% sim_init

delta_t1 = 0;
delta_t2 = 0;
amp_DeltaT1 = 50;
amp_DeltaT2 = 0;

sim ctrl_2x2 %Silumación no lineal

h(1) = subplot(2,2,1);
plot(t_sim,x_sim(:,params.op))
hold on
plot(t1,y1*amp_DeltaT1+params.Teq(params.op),'r')
plot(tlin,ylin(:,1,1)*amp_DeltaT1+params.Teq(params.op),'g')

h(3) = subplot(2,2,3); 
plot(t_sim,x_sim(:,params.op2))
hold on
plot(tlin,ylin(:,2,1)*amp_DeltaT1+params.Teq(params.op2),'g')

delta_t1 = 0;
delta_t2 = 0;
amp_DeltaT1 = 0;
amp_DeltaT2 = 50;
sim ctrl_2x2

h(2) = subplot(2,2,2); 
plot(t_sim,x_sim(:,params.op))
hold on
plot(tlin,ylin(:,1,2)*amp_DeltaT2+params.Teq(params.op),'g')

h(4) = subplot(2,2,4); 
plot(t_sim,x_sim(:,params.op2))
hold on
plot(t2,y2*amp_DeltaT2+params.Teq(params.op2),'r')
plot(tlin,ylin(:,2,2)*amp_DeltaT2+params.Teq(params.op2),'g')



%% Simulacion
close all
tend = 500;
delta_t1 = 150;
delta_t2 = 300;
amp_DeltaT1 = -1200;
amp_DeltaT2 = +500;
sim ctrl_2x2

figure
% salidas y acciones
subplot(4,1,1),plot(t_sim,x_sim(:,params.op))
subplot(4,1,2),plot(t_sim,p_sim)
subplot(4,1,3),plot(t_sim,x_sim(:,params.op2))
subplot(4,1,4),plot(t_sim,v_sim)

figure,plot_perfil
