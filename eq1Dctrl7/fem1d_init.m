function data = fem1d_init(X,params)
% Calculo de datos constantes para fem3d3

data.n=length(X);
data.nel=data.n-1;
data.Cconst=sparse(data.n,data.n);
data.Kconst=sparse(data.n,data.n);
data.fconst=zeros(data.n,1);
% data.fpoint=zeros(data.n,1);
data.nodes_point=[];
data.values_point=[];
[~,i]=min(abs(X-params.PotXt(0)));
el_point = i;
if (X(i)>params.PotXt(0) & i>1) el_point=i-1; end
data.X=X;
NN2=[2 1; 1 2]/6;
data.face_elements{1}=[1];
data.face_elements{2}=[length(X)-1];

for i=1:data.nel
   nodes=[i i+1];
    AA=[ones(2,1) X(nodes)'];
    lx=diff(X(nodes));
    abcd=AA\eye(2);
    B=[abcd(2,1) abcd(2,2)];
    data.Kek(:,:,i)=(B'*B)*lx*params.A;

    vu=[1 1]; % velocidad unitaria 1 o -1
    data.Kev(:,:,i)=(vu'*B)*lx*params.A*params.rho/2;
    data.Cel(i)=lx*params.A*params.rho/2; % para matriz C diagonal

    data.nodes(:,i)=nodes;
    data.cnodes(:,i)=int32(nodes)-1;
	
	data.Keh(:,:,i)=params.h*params.P*lx*NN2;
	data.Ker(:,i)=zeros(2,1); % para Kr (radiacion)
	data.feh(:,i)=params.h*params.P*params.Tinf*lx*[1; 1]/2;
    
    % entrada muy puntual
    if (params.point_input == 1 & i == el_point)
        N1=abcd(1,1)+abcd(2,1)*params.PotXt(0);
        N2=abcd(1,2)+abcd(2,2)*params.PotXt(0);
        data.nodes_point=nodes;
        data.values_point=[N1;N2];
    end
    
    data.Kconst(nodes,nodes) = data.Kconst(nodes,nodes) + data.Keh(:,:,i);
    data.fconst(nodes) = data.fconst(nodes) + data.feh(:,i);
end


end