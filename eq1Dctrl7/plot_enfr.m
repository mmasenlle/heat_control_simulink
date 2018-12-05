[~,xi] = min(abs(nodes - params.input_x));
dx=(x_sim(:,3:end) - x_sim(:,1:end-2))/(nodes(3)-nodes(1));
enf=dx./x_sim(:,2:end-1);

yy=enf(1:15:end,:);
tt=t_sim(1:15:end,:);

% close all
figure
g=plot(nodes(xi+1:end-1),zeros(size(nodes(xi+1:end-1))));
% m=min(min(yy)); M=max(max(yy))*1.2;
% hold on
% plot([params.input_x params.input_x],[m M],'m--');
% plot(nodes([params.op params.op]),[m M],'y--');
% plot(nodes([params.ep params.ep]),[m M],'g--');
% axis([0 params.Lx m M]);
ylim([-7 -1])
for i=1:size(yy,1)
    set(g,'ydata',yy(i,xi:end));
    title(['enfriamiento normalizado en t=',num2str(tt(i)),'s']);
    pause(.01);
end