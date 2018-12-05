function H = enthalpies_step(T,params)

% Entalp�as "reales": modelo sin "mushy region"

T0 = 0;
H0 = -7.98499888e+05; % Valores que establecen la referencia de entalp�a como 0 [J] en agua l�quida a 0 [�C]

% H = H0 + params.rho*params.cs*(T-T0); % solido
% H(find(T >= params.Tm)) = H0 + params.rho*(params.cs*(params.Tm-T0) + params.L + params.cl*(T(find(T >= params.Tm))-params.Tm));

H = H0 + params.rho*params.ccs*(T-T0); % solido
if T >= params.Tm
    H = H0 + params.rho*(params.ccs*(params.Tm-T0) + params.L + params.ccl*(T-params.Tm));
end

