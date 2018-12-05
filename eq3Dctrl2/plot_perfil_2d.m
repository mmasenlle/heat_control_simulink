%% Plot perfil 2d
TT=x_sim(1:1:end,:);
tt=t_sim(1:1:end,:);

cnodes=find(dtri.X(:,2) == params.input_y & dtri.X(:,3) == params.input_z);
g=plot(dtri.X(cnodes,1),TT(1,cnodes));
hold on
m=min(min(x_sim)); M=max(max(x_sim));
plot([params.input_x params.input_x],[m M],'m--');
plot([dtri.X(params.op,1) dtri.X(params.op,1)],[m M],'y--');
plot([dtri.X(params.op2,1) dtri.X(params.op2,1)],[m M],'c--');
for i=2:size(TT,1)
    set(g,'ydata',TT(i,cnodes));
    title(['perfil en t=',num2str(tt(i)),'s']);
    pause(.05);
end