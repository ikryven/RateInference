function [T Y] = integrate_macromodel( x0, t_end, coefs, params )

opt = odeset( 'AbsTol', 1e-25, 'RelTol', 1e-15 );

[ T Y ] = ode15s( @(t,x)  RHS_exp_poly( t, x, params.nu, params.S, coefs, params.C, x0(1), params.K ),  [0, t_end ], x0 ,opt);

end

function y = RHS_exp_poly( t, x, nu, S, coefs, F, x0, K );

    n = size( S, 1 );
    c = 1-x(1)/x0;
    x_nu = [];
    for j = 1:2
        x_nu( j, 1:size( x, 2 ) ) = prod( x.^( nu( j,: )' ) );
        k( j ) = F *exp( - polyval( coefs(j,:), c ) )* K{j};
    end
    

    for i = 1:length( x );
        y(i) = 0;
        for j = 1:n
             y(i) =  y(i) + k(j) * S( j, i ) * x_nu( j ); 
        end;
    end;
    y(2) = 0;
    y = y( : );

end
