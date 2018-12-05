function [sys,x0,str,ts,simStateCompliance] = ss1D_lin_2input(t,x,u,flag)
%
% Modelo para obtencion de linealizaciones potencia y velocidad
%

global iters ttimes params data;

switch flag,
    
    case 0, % Initialization %
        tic
        iters=0;
        ttimes=zeros(4,1);
        n=length(params.Teq);
        
        sizes = simsizes;
        
        sizes.NumContStates  = n;
        sizes.NumDiscStates  = 0;
        sizes.NumOutputs     = 2; % [T(p1),T(p2)]
        sizes.NumInputs      = 2; % [power,vel]
        sizes.DirFeedthrough = 0;
        sizes.NumSampleTimes = 1; % ts=0, offset=0
        sys = simsizes(sizes);

        x0  = params.Teq; %Estado inicial
        str = []; %Siempre es as�
        ts  = [0 0]; %Sist continuo
        simStateCompliance = 'UnknownSimState'; %Valor por defecto
        
        params.Pot=@(t) (params.Pot_nom);
        params.v= @(t) (params.v_nom); % m/s
        
    case 1, 
         params.Pot=@(t) (u(1));
         params.v=@(t) (u(2));
        
        sys = derTemp3_run(params,data,x,t);
        
    case 3, % Outputs %
        sys = [x(params.op); x(params.op2)];
    case 2, sys=[]; % Update: No hay estados discretos% 
        
    case 4, sys=[]; % GetTimeOfNextVarHit %
    case 9, sys=[]; % Terminate %

    otherwise, % Unexpected flags %
        DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
        
end
