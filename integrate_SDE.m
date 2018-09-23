function Y=integrate_SDE( tsim, coefs, X0, N_trajectories, params )
% X0 initial conditions

Y={};
Y{1} = [];
Y{2} = [];

for trajectory =1:N_trajectories
    disp( ['Trajectory # ' num2str(trajectory)] ); 
    x  = X0; 
    z  = [];
    yc = [];

    for i = 1:length( tsim ) - 1

        for j=1:size(X0)
            c = 1-x(1,i)/x(1,1);
            lambda = exp( - polyval( coefs(j,:), c ) );
             dt = tsim( i+1 ) - tsim( i );
           z( j, i )  = poissrnd(  lambda * prod( x( :, i ).^( params.nu( j,: )' ) ) * dt, 1, 1 );

        end;

        x(  :, i + 1 ) =  x( :, i ) + params.S * z( :, i );       % update
             
    end;
    
    Y{1}( end+1, 1:length( x ) ) = x( 1, : );
    Y{2}( end+1, 1:length( x ) ) = x( 2, : );
end


