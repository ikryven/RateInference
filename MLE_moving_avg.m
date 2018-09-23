function [ lambda var_lambda ] = MLE_moving_avg( X, Z, t, nu, s );

[n L] = size(Z);
L=L+1;
taus = diff( t ) ;

% x^nu
x_nu = [];
for j = 1:n
    x_nu( j, 1:size( X, 2 ) ) = prod( X.^( nu( j,: )' ) );
end


lambda = [];
var_lambda = [];
for j = 1:n
    % effective rate as a time series
    lambda( j, 1:L-1 ) = movmean( Z(j,:), s ) ./ movmean( x_nu( j, 1:L-1 ) .* taus, s );
    % varaince 
    var_lambda(j, 1:L-1 ) = lambda(j,:).^2    ./  movmean( Z(j,:), s ) / (2*s+1);
end;



