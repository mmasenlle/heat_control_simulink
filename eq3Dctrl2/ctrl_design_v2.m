% script para diseÃ±ar control de potencia
%% Linealización control de potencia
clear all; close all;
% para modificar parametros fisicos editar sim_init
tend=10;
sim_init;

% linealizaciones
model = 'lin_2in_3d';
opspec = operspec(model)
opspec.Inputs(1).u = params.Pot_nom;
opspec.Inputs(2).u = params.v_nom;
opspec.Inputs(1).Known = true;
opspec.Inputs(2).Known = true;
[op,opreport] = findop(model,opspec);
linsys = linearize(model,op);

% figure,step(linsys(1,1))
% figure,bode(linsys),grid
% figure,pzmap(linsys(1,1)), axis equal

%% Validación del modelo lineal
close all
figure
tend = 20; 
[ylin,tlin] = step(linsys,tend);

sim_init
t_step = 0;
delta_k = 0; delta_cp = 0; delta_pot = 100; delta_vel = 0;
sim open_loop_3d;

h(1) = subplot(2,2,1); 
plot(t_sim,x_sim(:,params.op),tlin,ylin(:,1,1)*delta_pot+x_sim(1,params.op))
h(3) = subplot(2,2,3); 
plot(t_sim,x_sim(:,params.op2),tlin,ylin(:,2,1)*delta_pot+x_sim(1,params.op2))

sim_init
t_step = 0;
delta_k = 0; delta_cp = 0; delta_pot = 0; delta_vel = 0.00005;
sim open_loop_3d;

h(2) = subplot(2,2,2); 
plot(t_sim,x_sim(:,params.op),tlin,ylin(:,1,2)*delta_vel+x_sim(1,params.op))
h(4) = subplot(2,2,4); 
plot(t_sim,x_sim(:,params.op2),tlin,ylin(:,2,2)*delta_vel+x_sim(1,params.op2))

%% Análisis RGA
close all
w = logspace(-3,2,100);
Pf = freqresp(linsys,w);
P_frd = frd(Pf,w);
% figure; bode(linsys,P_frd)

for i=1:length(w)
    RGAw(:,:,i)=Pf(:,:,i).*inv(Pf(:,:,i)).';
    RGAnum(:,:,i) = sum(sum(abs(RGAw(:,:,i) - eye(2))));
end

figure
subplot(2,2,1); semilogx(w,abs(squeeze(RGAw(1,1,:))))
subplot(2,2,2); semilogx(w,abs(squeeze(RGAw(1,2,:))))
subplot(2,2,3); semilogx(w,abs(squeeze(RGAw(2,1,:))))
subplot(2,2,4); semilogx(w,abs(squeeze(RGAw(2,2,:))))

RGAw(:,:,1)

%% QFT SISO

w_des = [0.001 0.01 0.1 1 10 100];
phs=0:-1:-360;
w1=logspace(-3,4,1000);

% Lazo 1
Wstab1 = 0.5/(cos(pi*(180-40)/(2*180)));
B1_stab = sisobnds(1,w_des,Wstab1,P_frd(1,1),[],1,[],1,phs);
P110 = P_frd(1,1,1);
% G110 = tf([1 1], [1 0]); % Manu
% G110 = tf(1, [1 0])*tf([1/2.1^2 2*0.7/2.1 1],[1/8.5^2 2*0.98/8.5 1]);
G110 = tf(1, [1 0])*tf([1/2^2 2*1/2 1],[1/8^2 2*0.75/8 1]);
lpshape(w1,B1_stab,P110,G110,phs);

% Lazo 2
Wstab2 = 0.5/(cos(pi*(180-40)/(2*180)));
B2_stab = sisobnds(1,w_des,Wstab1,P_frd(2,2),[],1,[],1,phs);
P220 = P_frd(2,2,1);
% G220 = tf(1,[1 0]); % Manu
% G220 = tf(-2.3e-6,[1 0]);
G220 = tf(-2.3e-6,[1 0])*tf([1/0.8 1],[1/2.8 1]);
lpshape(w1,B2_stab,P220,G220,phs);

G = [G110 tf(0,1); tf(0,1) G220];
T = feedback(linsys*G,eye(2));
step(T)

%% Tracking
sim_init
close

% Controladores de Manu
% Gc_pot = 1*tf([1 1], [1 0]);
% Gc_vel = 0*tf(1,[1 0]);

% Controladores loop-shaping
Gc_pot = tf(G110);
Gc_vel = tf(G220);

% Discretización
Gc_pot_d=c2d(Gc_pot,params.th,'tustin');
Gc_vel_d=c2d(Gc_vel,params.th,'tustin');

% Simulación
tend=20;
enable_ctrl_pot = 1;
enable_ctrl_vel = 1;
t_step_T1 = 0; amp_DeltaT1 = 200;  
t_step_T2 = 10; amp_DeltaT2 = 200;
t_step_k = 0; delta_k = 0;
t_step_cp = 0; delta_cp = 0;

sim ctrl_2x2

figure; set(gcf,'position',[680 49 560 948])
h(1) = subplot(4,1,1); plot(t_sim,x_sim(:,params.op)),ylabel('T1'),grid
h(2) = subplot(4,1,2); plot(t_sim,q_sim),ylabel('Pot'),grid
h(3) = subplot(4,1,3); plot(t_sim,x_sim(:,params.op2)),ylabel('T2'),grid
h(4) = subplot(4,1,4); plot(t_sim,qv_sim),ylabel('vel'),grid
linkaxes(h,'x')

[ylin,tlin] = lsim(T,[ref_T1 ref_T2],t_sim);
subplot(4,1,1); hold on; plot(tlin,ylin(:,1)+x_sim(1,params.op),'g')
subplot(4,1,3); hold on; plot(tlin,ylin(:,2)+x_sim(1,params.op2),'g')

% plot_0_end_x
% plot_0_end_y
% plot_perfil_2d
% plot_perfil
% plot_triple

%%
sys2x2_3_0_0_3 = linsys;
save sys2x2_3_0_0_3 sys2x2_3_0_0_3


%% QFT MIMO
for i=1:length(w)
    Pf_inv(:,:,i) = inv(Pf(:,:,i));
end
Psom_frd = frd(Pf_inv,w);

w_des = [0.001 0.01 0.1 1 10 100];
phs=0:-1:-360;

Wstab1 = 0.5/(cos(pi*(180-40)/(2*180)));
Wstab2 = 0.5/(cos(pi*(180-40)/(2*180)));
B1_stab = sisobnds(1,w_des,Wstab1,1/Psom_frd(1,1),[],1,[],1,phs);
B2_stab = sisobnds(1,w_des,Wstab2,1/Psom_frd(2,2),[],1,[],1,phs);

% Nominal plants
P110 = 1/Psom_frd(1,1,1);
P220 = 1/Psom_frd(2,2,1);

w1=logspace(-3,4,1000);

% G110 = tf([1 1], [1 0]); % Manu
G110 = tf(1,[1 0])*tf([1/2 1], [1/9 1]);
% G110 = tf(1, [1 0])*tf([1/2.1^2 2*0.7/2.1 1],[1/8.5^2 2*0.98/8.5 1]);
lpshape(w1,B1_stab,P110,G110,phs);

% G220 = tf(1,[1 0]); % Manu
% G220 = tf(-1,[1 0])*tf([1/650 1],[1/3200 1]); 
G220 = tf(-2.3e-6,[1 0]);
lpshape(w1,B2_stab,P220,G220,phs);

G = [G110 tf(0,1); tf(0,1) G220];
T = feedback(linsys*G,eye(2));
step(T)

%% %%%%% CÓDIGO DE MANU %%%%%

%% diseño controlador potencia
close all
%lpshape(linsys(1,1))
Gc=1*tf([1 1], [1 0]);
Tcl=feedback(Gc*linsys(1,1), 1);
figure,step(Tcl)
Tcl2=feedback(G110*linsys(1,1), 1);
figure,step(Tcl,Tcl2)

figure,bode(linsys(1,1),Gc),grid
legend('G','Gc');

% % Diseño para estabilidad
% w=1; % Como no hay incertidumbre, basta una frecuencia (cualquiera)
% W1=0.5/cos(130*pi/180/2); %MF=50º
% phs=-360:1:0;
% B_stab=sisobnds(1,w,W1,G,0,1,[],1,phs)
% plotbnds(B_stab,[],phs);
% title('Bounds de estabilidad');
% w1=logspace(-3,2,1000);
% lpshape(w1,B_stab,G,Gc,phs)


%% Validación del control mediante escalón
% clear all; close all; sim_init
Gc_pot=Gc;
Gc_pot_d=c2d(Gc_pot,params.th,'tustin');
% control de velocidad independiente
Gc_vel_d=c2d(tf(1,1),params.th,'tustin');
tend=30;
enable_ctrl_pot = 1;
enable_ctrl_vel = 0;
t_step_T1 = 5;
amp_DeltaT1 = 100;
t_step_T2 = 0;
amp_DeltaT2 = 0;
t_step_k = 0;
delta_k = 0;
t_step_cp = 0;
delta_cp = 0;

sim ctrl_2x2_invv

% plot_acciones
% plot_0_end_x
% plot_0_end_y
% plot_perfil_2d
% plot_perfil
% plot_triple

%% Simulación frente a cambios en velocidad
sim_init
% clear all; close all; sim_init
Gc_pot=Gc;
Gc_pot_d=c2d(Gc_pot,params.th,'tustin');
% control de velocidad independiente
Gc_vel_d=c2d(tf(1,1),params.th,'tustin');
tend=20;
enable_ctrl_pot = 1;
enable_ctrl_vel = 0;
t_step_T1 = 0;
amp_DeltaT1 = 0;
t_step_T2 = 0;
amp_DeltaT2 = 0;
t_step_k = 0;
delta_k = 0;
t_step_cp = 0;
delta_cp = 0;
% Ensayo en lazo abierto: cambiamos la velocidad
vold = params.v_nom;
params.v_nom = .001;

sim ctrl_2x2_invv
params.v_nom = vold;

% plot_acciones
% plot_0_end_x
% plot_0_end_y
% plot_perfil_2d
% plot_perfil
% plot_triple

%% Simulación frente a cambios en material
sim_init
Gc_pot=Gc;
Gc_pot_d=c2d(Gc_pot,params.th,'tustin');
% control de velocidad independiente
Gc_vel_d=c2d(tf(1,1),params.th,'tustin');
tend=25;
enable_ctrl_pot = 1;
enable_ctrl_vel = 0;
t_step_T1 = 0;
amp_DeltaT1 = 0;
t_step_T2 = 0;
amp_DeltaT2 = 0;
t_step_k = 0;
delta_k = 12;
t_step_cp = 0;
delta_cp = 0;

sim ctrl_2x2_invv

% plot_acciones
% plot_0_end_x
% plot_0_end_y
% plot_perfil_2d
% plot_perfil
% plot_triple



%% Control acoplado
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Linealización control de velocidad acoplado
clear all; close all;
% para modificar parametros fisicos editar sim_init
tend=10;
sim_init;

T1_ref = params.Teq(params.op);
Gc_pot=1*tf([1 1], [1 0]);


% linealizaciones
model = 'lin_v_ctrl_3d_invv';
opspec = operspec(model)
opspec.Inputs(1).u = 1/params.v_nom;
opspec.Inputs(1).Known = true;
[op,opreport] = findop(model,opspec);
linsys = linearize(model,op);

G3=linsys;
% figure, bode(G3), grid

%% diseño controlador de velocidad acoplado
close all
%lpshape(G3)
Gc=1*tf(1,[1 0]);

Tcl=feedback(Gc*G3, 1);
figure,bode(G3,Gc,Tcl),grid
figure,bode(G3*Gc),grid
legend('G3','Gc','Tcl');
%figure,step(G2,Tcl)
figure,step(Tcl),grid

%% Tracking
sim_init

% Controladores de Manu
Gc_pot = 1*tf([1 1], [1 0]);
Gc_vel = 1*tf(1,[1 0]);

% Controladores loop-shaping
Gc_pot = tf(G110);
% Gc_vel = tf(0,1);
Gc_vel = tf(G220);

% Discretización
Gc_pot_d=c2d(Gc_pot,params.th,'tustin');
Gc_vel_d=c2d(Gc_vel,params.th,'tustin');


% Simulación
tend=50;
enable_ctrl_pot = 1;
enable_ctrl_vel = 1;
t_step_T1 = 2; amp_DeltaT1 = 400;
t_step_T2 = 20; amp_DeltaT2 = 200;
t_step_k = 0; delta_k = 0;
t_step_cp = 0; delta_cp = 0;

sim ctrl_2x2

plot_acciones
% plot_0_end_x
% plot_0_end_y
% plot_perfil_2d
% plot_perfil
% plot_triple


%% Simulación frente a cambios en material
sim_init
Gc_pot=1*tf([1 1], [1 0]);
Gc_pot_d=c2d(Gc_pot,params.th,'tustin');
% control de velocidad independiente
Gc_vel=1*tf(1,[1 0]);
Gc_vel_d=c2d(Gc_vel,params.th,'tustin');
tend=80;
enable_ctrl_pot = 1;
enable_ctrl_vel = 1;
t_step_T1 = 0;
amp_DeltaT1 = 0;
t_step_T2 = 0;
amp_DeltaT2 = 0;
t_step_k = 5;
delta_k = 12;
t_step_cp = 0;
delta_cp = 0;

sim ctrl_2x2_invv

% plot_acciones
% plot_0_end_x
% plot_0_end_y
% plot_perfil_2d
% plot_perfil
% plot_triple

