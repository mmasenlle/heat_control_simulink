% clear all; close all;
global params nodes data
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
params.h=100000; % coeficiente de conveccion
params.a=0; % coeficiente de radiacion
params.Ta=0; % Temperatura inicial
params.Tinf=0; % Temperatura ambiente
params.lumped=1; % matriz de capacitancias diagonal
% dejar A=1 para que coincida con la funcion de transfereancia
params.A=1; % Area de la seccion de la barra unidimensional, sirver para escalar
params.Area=1; % Area de la seccion para la conveccion
params.P=2*pi*(sqrt(params.Area/pi)); % Perimetro de la barra unidimensional
% potencia nominal 
params.Pot_nom = 2e7; % (w)
% velocidad nominal
params.v_nom = .002; % (m/s)

% longitud de la pieza unidimensional
params.Lx=.4; % (m)
params.nx=201; % numero de nodos

% periodo para simulink
params.th = 0.05;

% punto de aplicacion de la fuente
params.input_x = params.Lx/8;
% punto donde se mide temperatura
params.output_x = params.input_x + (.01);
% punto donde se mide temperatura 2
params.output_x2 = params.input_x + (.11);
% punto donde se mide enfriamiento
params.enfr_x = params.input_x + (.005);


% valores iniciales calculados
nodes=linspace(0,params.Lx,params.nx);
% plot(nodes,ones(size(nodes)),'o');

params.Pot=@(t) (params.Pot_nom);
params.PotXt = @(t) (params.input_x);
params.PotWide= @(t) (0.002); % Anchura de la fuente de entrada (sigma de gauss)
params.point_input = 1; % indica fuente puntual
[~,params.op] = min(abs(nodes - params.output_x));
[~,params.op2] = min(abs(nodes - params.output_x2));
[~,params.ep] = min(abs(nodes - params.enfr_x));
params.v= @(t) (params.v_nom); % m/s
params.c_t = @(t) (params.cs);
params.k_t = @(t) (params.ks);
params.T0=ones(length(nodes),1)*params.Ta; % estado

% Nodos de contorno
params.faces{1}=[1];
params.faces{2}=[length(nodes)];

% propiedades en formato vector para cada nodo
params.ccs = ones(size(params.T0))*params.cs;
params.ccl = ones(size(params.T0))*params.cl;
params.kks = ones(size(params.T0))*params.ks;
params.kkl = ones(size(params.T0))*params.kl;

params.ind_mat_input = 1;

params.ceff_centro = (enthalpies_step(params.Tm+params.Teps,params) - enthalpies_step(params.Tm-params.Teps,params)) / (2*params.Teps) / params.rho;
params.Lmov1 = params.L/(params.Tm*params.Teps*2);
params.Lmov2 = params.Lmov1*(params.Tm - params.Teps);

params.Ts=0; % temp fija de contorno
params.ndirich=[params.faces{1}];
params.vdirich=ones(size(params.ndirich))*params.Ts;
params.T0(params.ndirich)=params.vdirich;

% calculo de datos geometricos y constantes
data=fem1d_init(nodes,params);

% busqueda del estado de equilibrio
[~,K,f] = fem1d_run(params,data,params.T0,0);
% dirichlet manteniendo K simetrica
bcwt = trace(K)/params.nx;
f = f - K(:,params.ndirich)*params.vdirich;
K(:,params.ndirich) = 0; K(params.ndirich,:) = 0;
K(params.ndirich,params.ndirich) = bcwt*speye(length(params.ndirich));
f(params.ndirich) = bcwt*params.vdirich;
% estado estacionario
params.Teq = K\f;

% Busqueda de equilibrio cuando hay cambio de fase
if (max(params.Teq) >= params.Tm & params.L > 0)
    err=1e30; alfa=0.9;
    for i=1:1000
        [~,K,f] = fem1d_run(params,data,params.Teq,0);
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

disp(sprintf('Init data en %f segundos', toc));

