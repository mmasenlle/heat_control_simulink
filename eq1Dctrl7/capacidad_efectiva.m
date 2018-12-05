function [cpe1 cpe2] = capacidad_efectiva(T,params)

% Capacidad efectiva calculada por el método de la sustitución directa

cpe1 = params.ccs;
cpe2 = params.ccs;
ind_med = find(T > (params.Tm-params.Teps) & T < (params.Tm+params.Teps));
ind_liq = find(T >= (params.Tm+params.Teps));
cpe1(ind_med) = params.ceff_centro(ind_med);
cpe1(ind_liq) = params.ccl(ind_liq);
cpe2(ind_med) = params.ccl(ind_med) + T(ind_med)*params.Lmov1 - params.Lmov2;
cpe2(ind_liq) = params.ccl(ind_liq) + params.L./T(ind_liq);