%demo
%% 
% Load species counts    -> X,
% reaction firing counts -> Z, 
% timesteps              -> t, 
% simultioan temperature -> T:
load MD_simulations/MD_T300K1.mat
params = simulation_parameters(T);


% output for reaction #
j=1;

%% = = = = = =  rate coefficient inferrence, time dependant
order = 3; 
%normalised time running from 0 to 1 for numerical stability
ct = t( 1: end-1 ) / max( t( 1: end-1 )  ); 
% all estimators
[ L_MLE1 VL_MLE1 ] = MLE_constant( X, Z, t, params.nu );
[ L_MLE2 VL_MLE2 ] = MLE_moving_avg( X, Z, t, params.nu, 1 );
[ alpha_t  V iH ]    = MLE_exp_poly( X, Z, t, params.nu, ct, order );




% plotting
figure(1)
cla
plot( t( [ 1, end ] ), params.C * L_MLE1( j, : ) * [ 1 1 ],  'Color', [ 0.2 0.8 0.2 ], 'LineWidth', 4 );
hold on
plot( t( [ 1 end ] ),( params.C *  L_MLE1( j, : ) - params.sigma * sqrt( params.C^2*VL_MLE1( j, : )  )) * [ 1 1 ], 'r--' )
plot( t( [ 1 end ] ), (params.C *  L_MLE1( j, : ) + params.sigma * sqrt( params.C^2*VL_MLE1( j, : )  )) * [ 1 1 ], 'r--' )
plot( t( 1:end - 1 ),  params.C  *   L_MLE2( j, : ), 'Color', [0.01 0.01 0.2] );
plot_var( t(1:end-1),  params.C  * exp( - polyval( alpha_t( j, : ), ct ) - params.sigma * sqrt( V(j,:)  ) ),  params.C* exp( - polyval( alpha_t(j,:), ct ) + params.sigma * sqrt( V(j,:)  ) ),'r' );
plot( t(1:end-1), params.C * exp( - polyval( alpha_t(j,:), ct ) ), 'r-', 'LineWidth', 2 );
xlabel('Time, $t$ [s]','Interpreter','Latex')
ylabel('Rate pre-factor, $A$ $[\frac{\mathrm{L}}{\mathrm{mol}\; \mathrm{s}}]$','Interpreter','Latex')
title([num2str(T) 'K'])
ylim([0 max(params.C * exp( - polyval( alpha_t( j, : ), ct )) + sqrt( V(j,:)  )) ] );




%% = = = = = =  rate coefficient inferrence, conversion dependant
order = 4; 
% conversion
cnv = 1 - X(1,:) / X(1,1);
[ lambda   ] = MLE_moving_avg( X, Z, t, params.nu, 1 );
[ alpha_chi,  V ] = MLE_exp_poly( X, Z, t, params.nu, cnv(1:end-1), order );


%plotting
figure(2)
cla
c = linspace(0,1,100)
stairs( cnv( 1:end - 1 ), params.C *   L_MLE2( j, : ), 'Color', [0.01 0.01 0.2] );
hold on
plot_var( cnv(1:end-1), params.C * exp( - polyval( alpha_chi( j, : ), cnv(1:end-1) ) - params.sigma * sqrt( V(j,:)  ) ),  params.C * exp( - polyval( alpha_chi(j,:), cnv(1:end-1) ) + params.sigma * sqrt( V(j,:)  ) ),'r' );
plot( c, params.C * exp( - polyval( alpha_chi(j,:), c ) ), '-', 'LineWidth', 2,'Color','r' );
xlabel('Conversion, $\chi$','Interpreter','Latex')
ylabel('Rate pre-factor, $A$ $[\frac{\mathrm{L}}{\mathrm{mol}\; \mathrm{s}}]$','Interpreter','Latex')
ylim([ 0 params.C * max( L_MLE2( j, : ) ) ] );
title([num2str(T) 'K'])
title([num2str(T) 'K, $' num2str(order)  '^\mathrm{th}$ order'],'Interpreter','Latex')
xlim( [0 1] )
set(gca,'yscale','lin')
set(gca,'xscale','lin')



%% = = = = = = Integration with SDE that reproduces Molecular Dynamics

%timesteps for SDE
tsim_sde = logspace( log10( t(2) ), log10( t(end) * 4e4 ), 300 );
N_trajectories = 100;
Y_SDE = integrate_SDE( tsim_sde, alpha_chi,  X( :, 1 ), N_trajectories, params );

%ploting
figure(3)
cla
plot( tsim_sde, Y_SDE{ 1 }', 'Color',[ 0.0 0 0.7 0.08 ], 'LineWidth',1 );
hold on
plot( t,    X( 1, : ), 'o', 'LineWidth', 5, 'Color', [ 0 1 1 ] );
plot( tsim_sde, Y_SDE{ 2 }', 'Color', [ 0.8 0 0.0 0.08], 'LineWidth',1 );
plot( t,    X( 2, : ), '+', 'LineWidth', 5, 'Color', [ 1 1 0 ] );

ylabel('\#V, \#R','Interpreter','Latex')
ylim([1 1e4])
xlim([1e-11 tsim_sde(end)])
xlabel('Time, $t$ [s]','Interpreter','Latex')
set(gca,'yscale','log')
set(gca,'xscale','log')


%% = = = = = = Integration with ODE that reproduces Molecular Dynamics

t_end = t(end) * 1e4; 
[t_sim_ode Y_ODE] = integrate_ODE_as_MD( t_end, alpha_chi, X(:,1), params );

%plotting
figure( 3 )
hold on

loglog( t_sim_ode, Y_ODE( :, 1 ) * params.Volume * params.N_a , '-k', 'Linewidth', 3 );
loglog( t_sim_ode, Y_ODE( :, 2 ) * params.Volume * params.N_a , '-k', 'Linewidth', 3 );

set( gca, 'yscale', 'log' )
set( gca, 'xscale', 'log' )
xlim( [ 1e-11 t_end ] )


%%
figure(4)
x0    =  X( :, 1 ) / params.C; %counts -> mol
x0(2) = 1e-3;

t_end = 100;
[T_macro Y_macro ] = integrate_macromodel( x0, t_end, alpha_chi, params );
loglog( T_macro, Y_macro( :, 1 )  , '-k', 'Linewidth', 3 );

set( gca, 'yscale', 'lin' );
set( gca, 'xscale', 'log' );
xlim( [ 1e-7 1e2 ] );

xlabel( 'Time, $t, [s]$', 'Interpreter', 'Latex' )
ylabel( 'Concentration, $c(t)$, $[\frac{\mathrm{mol} }{\mathrm{L}}]$','Interpreter','Latex' )


%% build table of rate constsnts
% Constant Rates in  [L/mol/s]
table_inferred_constatns()



