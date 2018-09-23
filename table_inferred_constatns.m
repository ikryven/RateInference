
batch{9} = { 'MD_simulations/MD_T600K1.mat', 'MD_simulations/MD_T600K2.mat', 'MD_simulations/MD_T600K3.mat', 'MD_simulations/MD_T600K4.mat'  };
batch{8} = { 'MD_simulations/MD_T550K1.mat', 'MD_simulations/MD_T550K2.mat', 'MD_simulations/MD_T550K3.mat'                   };
batch{7} = { 'MD_simulations/MD_T500K1.mat', 'MD_simulations/MD_T500K2.mat', 'MD_simulations/MD_T500K3.mat', 'MD_simulations/MD_T500K4.mat'  };
batch{6} = { 'MD_simulations/MD_T450K1.mat', 'MD_simulations/MD_T450K2.mat', 'MD_simulations/MD_T450K3.mat'                   };
batch{5} = { 'MD_simulations/MD_T400K1.mat', 'MD_simulations/MD_T400K2.mat', 'MD_simulations/MD_T400K3.mat', 'MD_simulations/MD_T400K4.mat'  };
batch{4} = { 'MD_simulations/MD_T350K1.mat', 'MD_simulations/MD_T350K2.mat', 'MD_simulations/MD_T350K3.mat'                   };
batch{3} = { 'MD_simulations/MD_T300K1.mat', 'MD_simulations/MD_T300K2.mat', 'MD_simulations/MD_T300K3.mat'                   };
batch{2} = { 'MD_simulations/MD_T250K1.mat', 'MD_simulations/MD_T250K2.mat'                                    };
batch{1} = { 'MD_simulations/MD_T200K2.mat', 'MD_simulations/MD_T200K3.mat'                                    };

TEMPS = [ 200 250 300 350 400 450 500 550 600 ]';

for id = 1:length( TEMPS )

    KK = [];
    VV = [];

    for fn = 1:length( batch{ id } );
        
        load( batch{id}{ fn } );
        
        params=simulation_parameters( T );
        
        [ L_MLE1 VL_MLE1 ] = MLE_constant( X, Z, t, params.nu );

        KK( end+1, : ) = params.C   *  [ L_MLE1(1)  * params.K{1}     L_MLE1(2)  * params.K{2}   ];
        VV( end+1, : ) = params.C^2 *  [ VL_MLE1(1) * params.K{1}^2  VL_MLE1(2)  * params.K{2}^2 ];

    end;

    k  = mean( KK, 1 ); 
    stdev = 2 * sqrt( sum( VV, 1 ) / size( VV, 1 )^2 );

    K_1(id)    = k(1);
    K_2(id)    = k(2);
    conf_1(id) = stdev(1);
    conf_2(id) = stdev(2);
    
    %LaTeX output:
    %disp([ num2str(TEMPS(id)) '  &  $' latex( vpa( k(1), 4 ) ) '$   & $\pm$ & $' latex( vpa( stdev(1), 4 ) ) '$    &    $' latex( vpa( k(2), 4) ) ' $ &  $\pm$ & $  ' latex( vpa( stdev(2), 4 ) ) '$ \\'])
    
end

K_1    = K_1(:);
K_2    = K_2(:);
conf_1 = conf_1(:);
conf_2 = conf_2(:);

TAB = table( TEMPS, K_1, conf_1, K_2, conf_2 )