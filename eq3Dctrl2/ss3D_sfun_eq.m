function [sys,x0,str,ts,simStateCompliance] = ss1D_sfun_eq(t,x,u,flag)
%
% Modelo an�logo el�ctrico en una dimensi�n de proceso de soldadura
%
% Entradas: u=[xd,cd,xa,ca]
%
% Estados: Las temperaturas de cada uno de los n nodos (n-1 intervalos iguales) en que se divide
%          la barra a soldar
%
% Salidas=Estados
%
% Es un sistema continuo en el tiempo
%

global iters ttimes params data f_old t_old;
% global mat1 mat2;

switch flag,
    
    case 0, % Initialization %
        tic
        iters=0;
        ttimes=zeros(4,1);
        n=length(params.Teq);
        
        sizes = simsizes;
        
        sizes.NumContStates  = 0;
        sizes.NumDiscStates  = n;
        sizes.NumOutputs     = n; % Todos los estados
        sizes.NumInputs      = 4; % [power,velocity,k(borde),cp(borde)]
        sizes.DirFeedthrough = 0;
        sizes.NumSampleTimes = 1; % ts=0, offset=0
        sys = simsizes(sizes);

        x0  = params.Teq; %Estado inicial
        str = []; %Siempre es as�
        ts  = [0 0]; %Sist continuo
        simStateCompliance = 'UnknownSimState'; %Valor por defecto
        
        params.Pot=@(t) (params.Pot_nom);
        params.v= @(t) (params.v_nom); % m/s
        params.k_t= @(t) (params.ks);
        params.c_t= @(t) (params.cs);
        
%         mat1=[];
%         mat2=[];
        t_old = 0;
        
        [~,~,f_old] = fem3d_run(params,data,x0,0);
        f_old(params.ndirich)=0; % dirichlet

    case 1, sys=[]; % Derivatives %
        
    case 3, % Outputs %
        sys = x;       
    case 2, % Update: No hay estados discretos% 
        params.Pot=@(t) (u(1));
        params.v = @(t) (u(2));
        params.k_t= @(t) (u(3));
        params.c_t= @(t) (u(4));
        
        theta = .5;
        Dt = params.th;
        T = x;

        iters=iters+1; tic1=tic;
        [C K f] = fem3d_run(params,data,T,t);
        lhs = Dt*(theta*f_old + (1-theta)*f);
        f(params.ndirich)=0; K(params.ndirich,:)=0; % dirichlet
        ttimes(1)=ttimes(1)+toc(tic1);tic2=tic;
        T=(C + theta*Dt*K)\((C - (1-theta)*Dt*K)*T + lhs);

        T(params.ndirich)=params.vdirich;
        
        ttimes(2)=ttimes(2)+toc(tic2);
        
        % movimiento de propiedades termicas
        mdt = params.Lx/(params.nx-1)/params.v(t);
        if (t_old + mdt < t)
            t_old = t;
            params.ccs(params.ind_c)=params.c_t(t);
            %params.ccl(params.ind_c)=params.c_t(t);
            params.kks(params.ind_c)=params.k_t(t);
            %params.kkl(params.ind_c)=params.k_t(t);
            
            params.ccs = params.MV*params.ccs;
            params.ccl = params.ccs; %params.ccl = MV*params.ccl;
            params.kks = params.MV*params.kks;
            params.kkl = params.kks; %params.kkl = MV*params.kkl;
        end
%         mat1=[mat1 params.ccs];
%         mat2=[mat2 params.kks];
        
        sys = T;
        f_old = f;
        
    case 4, sys=[]; % GetTimeOfNextVarHit %
    case 9, sys=[]; % Terminate %
        ttotal=toc;
        trest=ttotal-ttimes(1)-ttimes(2);
        str = sprintf('%d iteraciones en %f segundos', iters, ttotal);
        disp(str);
        str = sprintf('t1: %fs(%f%%), t2: %fs(%f%%), tr: %fs(%f%%)', ttimes(1),ttimes(1)*100/ttotal, ttimes(2),ttimes(2)*100/ttotal, trest,trest*100/ttotal);
        disp(str);

%anim_interfase;
    otherwise, % Unexpected flags %
        DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
        
end
