function derT = derTemp1_run(params,data,T,t)
global iters ttimes;
iters=iters+1;
tic1=tic;
[C K f] = fem1d_run(params,data,T,t);
f(params.ndirich)=0; K(params.ndirich,:)=0; % dirichlet
ttimes(1)=ttimes(1)+toc(tic1);
tic2=tic;
derT=C\(f-K*T);
% derT=sparse(1:data.n,1:data.n,C)\(f-K*T);
% [derT,flags]=bicg(C,(f-K*T));
% [derT,flags]=bicgstab(C,(f-K*T));
% [derT,flags]=bicgstabl(C,(f-K*T));
% [derT,flags]=cgs(C,(f-K*T));
% [derT,flags]=gmres(C,(f-K*T));
% [derT,flags]=lsqr(C,(f-K*T));
% [derT,flags]=qmr(C,(f-K*T));
% [derT,flags]=tfqmr(C,(f-K*T));
ttimes(2)=ttimes(2)+toc(tic2);
derT(params.ndirich)=0; % debería ser ya cero pero se fuerza por errores numéricos (FIXME)
end