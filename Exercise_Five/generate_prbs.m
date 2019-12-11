function [true_y] = generate_prbs(gamma, initial_state)
    A = [ 1 0 0 0 0 1 ;
          1 0 0 0 0 0 ;
          0 1 0 0 0 0 ;
          0 0 1 0 0 0 ;
          0 0 0 1 0 0 ;
          0 0 0 0 1 0 ] ;
    C = [0 0 0 0 0 2*gamma];
    x = []; 
    true_y = [] ;
    x(:,1) = initial_state' ;
    num_state_params = 6;
    max_period = (2^num_state_params-1);
    num_of_periods = 3;
    for i = 1:(num_of_periods*max_period) + 1
        z = A*x(:,i) ;
        x(:,i+1) = mod(z,2) ;
        true_y(i) = C*x(:,i) - gamma ; 
    end
end