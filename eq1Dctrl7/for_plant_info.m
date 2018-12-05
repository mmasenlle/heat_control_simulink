
desps=[15,20,25,30,40,50,60,80,100,120];

%desps=[30,60,100];

xi_mm=10;
xj_mm=15;

sim_init2
plant_info

figure
%subplot(1,2,1)
semilogx(w,rgas(:,1));
%subplot(1,2,2)
%semilogx(w,rgas(:,2));



rgas=[];
t_sim1=[];t_sim2=[];t_sim3=[];t_sim4=[];
x_sim1=[];x_sim2=[];x_sim3=[];x_sim4=[];
fr1=[];fr2=[];fr3=[];fr4=[];ww=[];
for i=1:length(desps)
    xi_mm=10;
    xj_mm=desps(i);
    
    sim_init2
    plant_info
end
figure,plot(desps,rgas(:,1),desps,rgas(:,2));
fname=sprintf('img/rgas_xi_%d.png', xi_mm);
saveas(gcf,fname);
close all

figure
h(1) = subplot(2,2,1); 
plot(t_sim1,x_sim1)
h(3) = subplot(2,2,3); 
plot(t_sim3,x_sim3)
h(2) = subplot(2,2,2); 
plot(t_sim2,x_sim2)
h(4) = subplot(2,2,4); 
plot(t_sim4,x_sim4)

fname=sprintf('img/step_xi_%d.png', xi_mm);
saveas(gcf,fname);
close all

figure
h(1) = subplot(2,2,1); 
semilogx(ww,20*log10(fr1))
h(3) = subplot(2,2,3); 
semilogx(ww,20*log10(fr3))
h(2) = subplot(2,2,2); 
semilogx(ww,20*log10(fr2))
h(4) = subplot(2,2,4); 
semilogx(ww,20*log10(fr4))

fname=sprintf('img/bode_xi_%d.png', xi_mm);
saveas(gcf,fname);
close all
