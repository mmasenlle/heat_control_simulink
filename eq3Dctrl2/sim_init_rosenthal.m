% clear all; close all;
global params dtri data
tic

% constantes termicas de un acero
% capacidad calorifica (solido y liquido)
params.cs = 460;    % J/kg/K
params.cl = 460; % J/kg/K
% conductividad (solido y liquido)
params.ks = 24.0; % J/m/s/K
params.kl = 24.0;  % J/m/s/K
% densidad, temperatura fusion y calor especifico
params.rho = 7925; % kg/m^3
params.Tm = 1783; % K
params.L = 0; %270000; % J/kg
% rango de temperatura durante la fusion
params.Teps = 30;
params.h=0; %1000; % coeficiente de conveccion
params.a=0; %0.3; % coeficiente de radiacion
params.Ta=0; % Temperatura inicial
params.Tinf=0; % Temperatura ambiente
% potencia nominal 
params.Pot_nom = 500; % (w)
% velocidad nominal
params.v_nom = .002; % (m/s)

% dimensiones de la pieza
params.Lx=.08; % (m)
params.Ly=.04; % (m)
params.Lz=.001; % (m)
params.nx=21; % numero de nodos en eje x
params.ny=11; % numero de nodos en eje y
params.nz=3; % numero de nodos en eje z

% periodo para simulink
params.th = 0.05;

% punto de aplicacion de la fuente
params.input_x = params.Lx/4;
params.input_y = params.Ly/2;
params.input_z = 0;
% punto donde se mide temperatura
params.output_x = params.input_x + (.003);
params.output_y = params.input_y;
params.output_z = params.input_z;
params.output_x2 = params.input_x;
params.output_y2 = params.input_y + (.003);
params.output_z2 = params.input_z;


% valores iniciales calculados
% mallado
X=linspace(0,params.Lx,params.nx);
Y=linspace(0,params.Ly,params.ny);
Z=linspace(0,params.Lz,params.nz);
[XX YY ZZ]=ndgrid(X,Y,Z);
n=params.nx*params.ny*params.nz;
XX=reshape(XX,[n 1]);YY=reshape(YY,[n 1]);ZZ=reshape(ZZ,[n 1]);
nodes=[XX YY ZZ]; clear XX YY ZZ; % plot3(nodes(:,1),nodes(:,2),nodes(:,3),'o'); axis equal
% triangulacion
dtri=DelaunayTri(nodes); % tetramesh(dtri); axis equal
nodes=dtri.X;

params.Pot=@(t) (params.Pot_nom);
% puntos de medicion
[~,xo] = min(abs(X - params.output_x)); % X(xo)
[~,yo] = min(abs(Y - params.output_y)); % Y(yo)
[~,zo] = min(abs(Z - params.output_z)); % Z(zo)
params.op = ((zo-1)*params.ny + (yo-1))*params.nx + xo; % dtri.X(params.op,:)
[~,xo] = min(abs(X - params.output_x2)); % X(xo)
[~,yo] = min(abs(Y - params.output_y2)); % Y(yo)
[~,zo] = min(abs(Z - params.output_z2)); % Z(zo)
params.op2 = ((zo-1)*params.ny + (yo-1))*params.nx + xo; % dtri.X(params.op2,:)

params.v= @(t) (params.v_nom); % m/s
params.c_t = @(t) (params.cs);
params.k_t = @(t) (params.ks);
params.T0=ones(n,1)*params.Ta; % estado

% Nodos de contorno
params.faces{1}=find(nodes(:,1)==0);
params.faces{2}=find(nodes(:,1)==params.Lx);
params.faces{3}=find(nodes(:,2)==0);
params.faces{4}=find(nodes(:,2)==params.Ly);
params.faces{5}=find(nodes(:,3)==0);
params.faces{6}=find(nodes(:,3)==params.Lz);
% face=1; close all; plot3(nodes(:,1),nodes(:,2),nodes(:,3),'o'); axis equal; hold on; plot3(nodes(params.faces{face},1),nodes(params.faces{face},2),nodes(params.faces{face},3),'ro'); clear face

% preparacion de matriz de movimiento de propiedades
mv1=sparse([1 zeros(1,params.nx-1); diag(ones(1,params.nx-1)) zeros(params.nx-1,1)]);
params.MV=mv1;
for i=2:(params.ny*params.nz)
    params.MV=blkdiag(params.MV,mv1);
end   % full(params.MV)
params.ind_c=1:params.nx:n; % norm(params.ind_c'-params.faces{1})
% close all; plot3(nodes(:,1),nodes(:,2),nodes(:,3),'o'); axis equal; hold on; plot3(nodes(params.ind_c,1),nodes(params.ind_c,2),nodes(params.ind_c,3),'ro');
% close all; x=zeros(n,1); x(params.ind_c)=1; x2=params.MV*x; x2(params.ind_c)=0; ind=find(x2); plot3(nodes(:,1),nodes(:,2),nodes(:,3),'o'); axis equal; hold on; plot3(nodes(ind,1),nodes(ind,2),nodes(ind,3),'ro'); clear x x2 ind

% limpieza de variables temporales
clear X Y Z nodes xo yo zo mv1 n;

% propiedades en formato vector para cada nodo
params.ccs = ones(size(params.T0))*params.cs;
params.ccl = ones(size(params.T0))*params.cl;
params.kks = ones(size(params.T0))*params.ks;
params.kkl = ones(size(params.T0))*params.kl;

params.ceff_centro = (enthalpies_step(params.Tm+params.Teps,params) - enthalpies_step(params.Tm-params.Teps,params)) / (2*params.Teps) / params.rho;
params.Lmov1 = params.L/(params.Tm*params.Teps*2);
params.Lmov2 = params.Lmov1*(params.Tm - params.Teps);

params.Ts=0; % temp fija de contorno
params.ndirich=[params.faces{1}];
params.vdirich=ones(size(params.ndirich))*params.Ts;
params.T0(params.ndirich)=params.vdirich;

% calculo de datos geometricos y constantes
data=fem3d_init(dtri,params);

% busqueda del estado de equilibrio
[~,K,f] = fem3d_run(params,data,params.T0,0);
% dirichlet manteniendo K simetrica
if (params.ndirich)
    bcwt = trace(K)/params.nx;
    f = f - K(:,params.ndirich)*params.vdirich;
    K(:,params.ndirich) = 0; K(params.ndirich,:) = 0;
    K(params.ndirich,params.ndirich) = bcwt*speye(length(params.ndirich));
    f(params.ndirich) = bcwt*params.vdirich;
end
% estado estacionario
params.Teq = K\f;

% Busqueda de equilibrio cuando hay cambio de fase
if (max(params.Teq) >= params.Tm & params.L > 0)
    err=1e30; alfa=0.9;
    for i=1:1000
        [~,K,f] = fem3d_run(params,data,params.Teq,0);
        f = f - K(:,params.ndirich)*params.vdirich;
        K(:,params.ndirich) = 0; K(params.ndirich,:) = 0;
        K(params.ndirich,params.ndirich) = bcwt*speye(length(params.ndirich));
        f(params.ndirich) = bcwt*params.vdirich;
        T = params.Teq*alfa + (1-alfa)*(K\f);
        er1=norm(T-params.Teq);
        if (er1 >= err) break; end
        params.Teq = T; err = er1;
    end
    disp(sprintf('Err: %g at #%d', err, i));
end
clear bcwt K f i;
disp(sprintf('Init data en %f segundos', toc));
% return

x_sim=[];x_sim(1,:)=params.Teq; t_sim=[0];
figure(1); plot_perfil;
figure(2); plot_triple;
x_sim %Temperaturas en equilibrio
dtri.X %coordenadas





