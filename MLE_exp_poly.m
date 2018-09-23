function [ alpha V iH EXITFLAG] = MLE_exp_poly( X, Z, t, nu, ct, s_max );
 s = 1; % for intitial guess;

[n L]    = size( Z );
L = L+1;
taus     = diff( t ) ;

alpha   = [];
x_nu = [];

for j = 1:n
    x_nu( j, 1:size( X, 2 ) ) = prod( X.^( nu( j,: )' ) );
end
x_nu( :, end ) = [];

%% time-series MLE needed for intitial guess
lambda = [];
var_lambda = [];
for j = 1:n
    
    % effective rate as a time series
    lambda( j, 1:L-1 ) = movmean( Z(j,:), s ) ./ movmean( x_nu( j, 1:L-1 ) .* taus, s );
    
    % varaince 
    var_lambda(j, 1:L-1 ) = lambda( j, : ).^2    ./  movmean( Z( j, : ), s ) / L;
    
end;


%% 
for j  = 1:n;

    id = ( lambda( j, : ) > 0 ) & ~isinf( lambda( j, : ) );

    for order=s_max:-1:0

        alpha0 = polyfit( ct( id ), -log( lambda( j, id ) ), order );
        if max(exp(-polyval( alpha0, ct )))<1e15
            break
        else
        %%%
        warning('Reducing the power of the initial guess')
        end

    end
    
    alpha0 = [zeros(1,s_max-order) alpha0];

    options = optimoptions( 'fsolve', 'Display', 'iter', 'StepTolerance', 1e-10, 'OptimalityTolerance', 1e-10, 'MaxFunctionEvaluations', 1e5, 'MaxIterations', 500, 'FunctionTolerance', 1e-10 );
    [ a, ~, EXITFLAG] = fsolve( @( alpha )  RHS( alpha, ct, x_nu( j, : ), taus, Z( j, : ) ), alpha0, options );

    alpha(j,1:s_max+1)=a;
    if alpha0(1)<0 & ct(end)>100
        warning('The leading coefficient is negatinve, are you using even polynomial power with large time?')
    end

end

%% confidence intervals
V = 0;
    for j = 1:n

        H=0;
        for k  = 0:s_max;
            for s  = 0:s_max;
                  H( k+1, s+1 )  = L * sum(  exp( - polyval( alpha( j, : ), ct ) ) .* ct.^s .* ct.^k .* x_nu( j, : ) .* taus  ); 
             end;
        end;
        iH{j} = inv(H);

        for l=1:length( ct )
            b = ct( l ) .^( s_max:-1:0 ); 
            
            V( j, l ) = b * iH{j} * b' ;
        end

    end



end

function y = RHS( alpha, ct, x_nu, taus, z );
    y=[];
    
    for s = 0:length(alpha)-1
        
        y( s + 1 ) =  sum( ( exp( - polyval( alpha, ct ) )  .* x_nu .* taus - z ) .* ct.^s ); 

    end; 
    
end


