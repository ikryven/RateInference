function [ lambda var_lambda ] = MLE_constant( X, Z, t, nu );


taus = [ diff( t )  0];

L = length(t);

% x^nu
x_nu = [];
for j = 1:2
    x_nu( j, 1:size( X, 2 ) ) = prod( X.^( nu( j,: )' ) );
end


% apparent rate constant:
lambda = mean( Z, 2 ) ./ mean( x_nu .* repmat( taus, 2, 1), 2 );

% varaince 
var_lambda = lambda.^2 ./ mean( Z, 2 ) / L;