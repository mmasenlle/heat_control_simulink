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

G=linsys(1,1);
% step(linsys(2,1))
% figure,step(G)
% figure,bode(G),grid
% figure,pzmap(G), axis equal

%% Validación del modelo lineal
sim_init
tend = 20;
t_step = 0;
delta_k = 0;
delta_cp = 0;
delta_pot = 50;
delta_vel = 0;
sim open_loop_3d;
figure;
[ylin,tlin] = step(G*delta_pot,tend);
plot(t_sim,x_sim(:,params.op),tlin,ylin+x_sim(1,params.op))

%% diseño controlador potencia
close all
%lpshape(G)
Gc=1*tf([2 1], [1 0]);
Tcl=feedback(Gc*G, 1);
figure,step(Tcl)
figure,bode(Gc*G),grid
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
tend=50;
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
Gc_pot=1*tf([2 1], [1 0]);


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
Gc=2*tf([1 1],[1 0]);

Tcl=feedback(Gc*G3, 1);
%figure,bode(G3,Gc,Tcl),grid
figure,bode(G3*Gc),grid
legend('G3','Gc','Tcl');
%figure,step(G2,Tcl)
figure,step(Tcl),grid

%% Tracking
sim_init
Gc_pot=1*tf([2 1], [1 0]);
Gc_pot_d=c2d(Gc_pot,params.th,'tustin');
% control de velocidad dependiente
Gc_vel=2*tf([1 1],[1 0]);
Gc_vel_d=c2d(Gc_vel,params.th,'tustin');
tend=50;
enable_ctrl_pot = 1;
enable_ctrl_vel = 1;
t_step_T1 = 2;
amp_DeltaT1 = 400;
t_step_T2 = 20;
amp_DeltaT2 = -100;
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


%% Simulación frente a cambios en material
sim_init
Gc_pot=1*tf([2 1], [1 0]);
Gc_pot_d=c2d(Gc_pot,params.th,'tustin');
% control de velocidad independiente
Gc_vel=2*tf([1 1],[1 0]);
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

