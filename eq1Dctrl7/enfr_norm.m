function y = enfr_norm(x,point)
% enfriamiento normalizado en punto point
global nodes
y = (x(point+1)-x(point-1))/(nodes(point+1)-nodes(point-1))/x(point);