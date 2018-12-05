%% Plot antes y después de la parte central
cnodes=find(dtri.X(:,2) == params.input_y & dtri.X(:,3) == params.input_z);
plot(dtri.X(cnodes,1),x_sim(1,cnodes),dtri.X(cnodes,1),x_sim(end,cnodes));
hold on
m=min(min(x_sim)); M=max(max(x_sim));
plot([params.input_x params.input_x],[m M],'m--');
plot([dtri.X(params.op,1) dtri.X(params.op,1)],[m M],'y--');
plot([dtri.X(params.op2,1) dtri.X(params.op2,1)],[m M],'c--');