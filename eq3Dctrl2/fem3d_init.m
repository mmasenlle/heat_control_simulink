function data = fem3d_init(dt,params)
% Calculo de datos constantes para fem2d

data.n=size(dt.X,1);
data.nel=size(dt.Triangulation,1);
data.Cconst=sparse(data.n,data.n);
data.Kconst=sparse(data.n,data.n);
data.fconst=zeros(data.n,1);
% data.fpoint=zeros(data.n,1);
data.nodes_point=[];
data.values_point=[];
el_point = pointLocation(dt, [params.input_x,params.input_y,params.input_z]);
data.dt=dt;

for i=1:data.nel
    nodes=dt.Triangulation(i,:);
    AA=[ones(4,1) dt.X(nodes,:)];
    V=det(AA)/6;
    abcd=AA\eye(4);
    B=[abcd(2,1) abcd(2,2) abcd(2,3) abcd(2,4)
       abcd(3,1) abcd(3,2) abcd(3,3) abcd(3,4)
       abcd(4,1) abcd(4,2) abcd(4,3) abcd(4,4)];
    data.Kek(:,:,i)=(B'*B)*V;
    vu=[1 1 1 1; 0 0 0 0; 0 0 0 0]; % velocidad unitaria en direccion x
    data.Kev(:,:,i)=(vu'*B)*V*params.rho/4;
    data.Cel(i)=V*params.rho/4; % para matriz C diagonal

    data.nodes(:,i)=nodes;
    data.cnodes(:,i)=int32(nodes)-1;
    
    % entrada muy puntual
    if (i == el_point)
        N1=abcd(1,1)+abcd(2,1)*params.input_x+abcd(3,1)*params.input_y+abcd(4,1)*params.input_z;
        N2=abcd(1,2)+abcd(2,2)*params.input_x+abcd(3,2)*params.input_y+abcd(4,2)*params.input_z;
        N3=abcd(1,3)+abcd(2,3)*params.input_x+abcd(3,3)*params.input_y+abcd(4,3)*params.input_z;
        N4=abcd(1,4)+abcd(2,4)*params.input_x+abcd(3,4)*params.input_y+abcd(4,4)*params.input_z;
        data.nodes_point=nodes;
        data.values_point=[N1;N2;N3;N4];
    end

% sin conveccion ni radiacion por el momento  
%    data.Kconst(nodes,nodes) = data.Kconst(nodes,nodes) + data.Keh(:,:,i);
%    data.fconst(nodes) = data.fconst(nodes) + data.feh(:,i);
end

end