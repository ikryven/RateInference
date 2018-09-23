function [T Y] = integrate_ODE_like_MD( t_end, coefs, x0, params )

x0    =  x0 / params.Volume / params.N_a; %convert integer counts to concentrations


opt = odeset( 'AbsTol', 1e-25 );
[ T Y ] = ode15s( @(t,x)  RHS_exp_poly( t, x, params.nu, params.S, coefs, params.Volume * params.N_a, x0(1) ),  [0, t_end ], x0 ,opt);


end

function y = RHS_exp_poly( t, x, nu, S, coefs, C, x0 );

    n = size( S, 1 );
    c = 1-x(1)/x0;
    x_nu = [];
    for j = 1:2
        x_nu( j, 1:size( x, 2 ) ) = prod( x.^( nu( j,: )' ) );
        k( j ) = C * exp( - polyval( coefs(j,:), c ) );
    end
    

    for i = 1:length( x );
        y(i) = 0;
        for j = 1:n
             y(i) =  y(i) + k(j) * S( j, i ) * x_nu( j ); 
        end;
    end;
    
    y = y( : );

end