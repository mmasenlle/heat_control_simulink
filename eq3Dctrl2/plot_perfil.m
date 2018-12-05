TT=x_sim(1:10:end,:);
tt=t_sim(1:10:end,:);

% close all
figure
mm1=dtri.X(params.faces{5},1); mm2=dtri.X(params.faces{5},2);
[mmg1,mmg2] = meshgrid(linspace(min(mm1),max(mm1),params.nx),linspace(min(mm2),max(mm2),params.ny));
mmgT = griddata(mm1,mm2,params.Teq(params.faces{5}),mmg1,mmg2);
g=surfc(mmg1,mmg2,mmgT); %view(0,90); axis equal
m = min(min(TT)); M = max(max(TT))*1.2;
axis([0 params.Lx 0 params.Lx m M]);
hold on
plot3([params.input_x params.input_x],[params.input_y params.input_y],[m M],'m--');
plot3([dtri.X(params.op,1) dtri.X(params.op,1)],[dtri.X(params.op,2) dtri.X(params.op,2)],[m M],'y--');
plot3([dtri.X(params.op2,1) dtri.X(params.op2,1)],[dtri.X(params.op2,2) dtri.X(params.op2,2)],[m M],'c--');
for i=1:size(TT,1)
    set(g,'zdata',griddata(mm1,mm2,TT(i,params.faces{5}),mmg1,mmg2));
    title(['perfil en t=',num2str(tt(i)),'s']);
    pause(.01);
end
return
