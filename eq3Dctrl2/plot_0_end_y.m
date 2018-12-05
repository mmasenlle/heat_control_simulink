%% Plot antes y después del corte transversal
tnodes=find(dtri.X(:,1) == dtri.X(params.op2,1) & dtri.X(:,3) == params.input_z);
plot(dtri.X(tnodes,2),x_sim(1,tnodes),dtri.X(tnodes,2),x_sim(end,tnodes));
hold on
m=min(min(x_sim)); M=max(max(x_sim));
plot([params.input_y params.input_y],[m M],'m--');
plot([dtri.X(params.op2,2) dtri.X(params.op2,2)],[m M],'y--');