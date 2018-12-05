TT=x_sim(1:50:end,:);
tt=t_sim(1:50:end,:);

% close all
figure
% set(gcf,'position',[300 200 600 700])

m1=dtri.X(params.faces{5},1);L1=max(m1);
m2=dtri.X(params.faces{5},2);L2=max(m2);
[mg1,mg2] = meshgrid(linspace(0,L1,50),linspace(0,L2,50));
mT = griddata(m1,m2,TT(1,params.faces{5}),mg1,mg2);
temp_levels = [0:200:1600 params.Tm-params.Teps params.Tm+params.Teps];

h(1)=subplot(3,1,1,'position',[0.1 0.7 0.8 0.25]);
g=surfc(mg1,mg2,mT);
axis([0 params.Lx 0 params.Ly 0 3000]);
view(0,0)
hold on;

h(2)=subplot(3,1,2,'position',[0.1 0.4 0.8 0.25]);
[~,g_sup] = contourf(mg1,mg2,mT,temp_levels);
caxis([0 2000]);
hold on;

cnodes=find(dtri.X(:,2) == params.input_y);
mm1=dtri.X(cnodes,1);Lm1=max(mm1);
mm2=dtri.X(cnodes,3);Lm2=max(mm2);
[mmg1,mmg2] = meshgrid(linspace(0,Lm1,50),linspace(0,Lm2,50));
mmt = griddata(mm1,mm2,TT(1,cnodes),mmg1,mmg2);
h(3) = subplot(3,1,3,'position',[0.1 0.1 0.8 0.25]);
[~,mg] = contourf(mmg1,-mmg2,mmt,temp_levels);
caxis([0 2000]);
hold on;

linkaxes(h,'x')

for i=2:length(tt)
    zdata=griddata(m1,m2,TT(i,params.faces{5}),mg1,mg2);
    %title(['t: ', num2str(tt(i)), ' seg']);
    set(g,'zdata',zdata);
    set(g_sup,'zdata',zdata);
    set(mg,'zdata',griddata(mm1,mm2,TT(i,cnodes),mmg1,mmg2));
    pause(0.1);
end
