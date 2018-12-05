% Rosenthal
alpha = params.ks/(params.rho*params.cs);
% tmax = alpha/params.v_nom/params.ks*params.Pot_nom;
m=sqrt((params.v_nom/(2*alpha))^2 + params.P*params.h/(params.A*params.ks));
tmax=params.Pot_nom/(2*params.A*params.ks*m);
x = nodes(1:end);
[~,xi] = min(abs(nodes + params.input_x - params.Lx));
x1 = fliplr(nodes(xi+1:end) + params.input_x - params.Lx); % distancia desde los puntos a la izquierda de la torcha hasta la torcha
x2 = fliplr(params.Lx - nodes(1:xi) - params.input_x); % distancia desde los puntos a la derecha de la torcha hasta la torcha

y = [tmax*exp((-m-params.v_nom/(2*alpha))*x1) tmax*exp((-m+params.v_nom/(2*alpha))*x2)];

% tmax = alpha/u0/params.ks*1e6;
% y = [tmax*exp(-params.v_nom/alpha*(x(end)/2 - x1)) tmax tmax*ones(size(x2))];

hold on;
plot(x,y,'r')