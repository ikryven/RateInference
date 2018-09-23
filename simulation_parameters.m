function params = simulation_parameters( T );


params.N_a     = 6.022e23;  % mol^-1
params.R       = 8.3144598; % J K^1 mol^1
params.m3_to_L = 1e3; %converts m3 to liters,  
params.sigma   = 2; % for confidence intervals

params.Volume = 7.5213e-25; % m^-3

%0.5 is the parameter used in MD ( probability that reaction fires reactants are within reaction radius ) 
params.P_react = 0.5;


% scaling constant in the pre-factor A(t)
params.C = params.Volume * params.N_a * params.m3_to_L / params.P_react;

% stoichiometry:
params.nu = [ 1 1; 0 2 ];
params.S  = [ -1   0  ;   0   -2 ];

%Energy for HDDA polymerisation
E_p{ 1 } = 30320.07 + 702; % RMGpy rules + the dip in Lennard Jones
E_p{ 2 } = 7971     + 702; % RMGpy rules + the dip in Lennard Jones

params.K{ 1 } = exp( - E_p{ 1 } / ( params.R * T ) ); % propagation
params.K{ 2 } = exp( - E_p{ 2 } / ( params.R * T ) ); % termination
params.T=T;