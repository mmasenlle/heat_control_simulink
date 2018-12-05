TT=x_sim(1:25:end,:);
tt=t_sim(1:25:end,:);

% close all
% figure
cla
g=plot(nodes,zeros(size(nodes)));
m=min(min(TT)); M=max(max(TT))*1.2;
hold on
plot([params.input_x params.input_x],[m M],'m--');
plot(nodes([params.op params.op]),[m M],'y--');
plot(nodes([params.op2 params.op2]),[m M],'c--');
%plot(nodes([params.ep params.ep]),[m M],'g--');
plot(nodes,x_sim(1,:),'k');
axis([0 params.Lx m M]);
for i=1:size(TT,1)
    set(g,'ydata',TT(i,:));
    title(['perfil en t=',num2str(tt(i)),'s']);
    pause(.01);
end