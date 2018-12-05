function [C K f] = fem3d_run(params,data,T,t)
% Creación de las matrices C K y el vector f dados
% los datos inicializados
global ttimes;

tic1=tic;
k=params.kks;
ind_liq=find(T > params.Tm);
k(ind_liq)=params.kkl(ind_liq);
[cpe1 cpe2]=capacidad_efectiva(T, params);
% plot(data.X,cpe1/10,data.X,cpe2);pause(.1);
[C K f] = fem3d_runc(data,params,T,t,cpe1,cpe2,k,params.v(t));
C=C+data.Cconst; K=K+data.Kconst; f=f+data.fconst;
f(data.nodes_point)=f(data.nodes_point)+(data.values_point*params.Pot(t));
return;

C=data.Cconst;
K=data.Kconst;
f=data.fconst;
f(data.nodes_point)=(data.values_point*params.Pot(t));
for i=1:data.nel
    C(data.nodes(:,i),data.nodes(:,i)) = C(data.nodes(:,i),data.nodes(:,i)) + ...
        data.Cel(i)*diag(cpe1(data.nodes(:,i)));
    K(data.nodes(:,i),data.nodes(:,i)) = K(data.nodes(:,i),data.nodes(:,i)) + ...
        data.Kek(:,:,i)*mean(k(data.nodes(:,i)))+data.Kev(:,:,i)*diag(cpe2(data.nodes(:,i)))*params.v(t);
end
end