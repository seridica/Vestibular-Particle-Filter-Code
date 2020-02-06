%%%
% File: ExplicitIntegrator.m
% Author: Calvin Kuo
% Date: 10-23-2018
% Notes: This code is for explicit integration related to particle filter
% project. Of note, this integrator can store values for later use.
% code originally written by PAF, CJD, NK, and JSB

function [tinteg, integ_states, store_states] = ExplicitIntegrator( diff_eq, tvec, states_init, store_init )
    
    % Setup outputs
    integ_states = zeros( length(tvec), length( states_init ) );
    store_states = zeros( length(tvec), length( store_init ) );

    integ_states(1,:) = states_init;
    store_states(1,:) = store_init;
    
    % Go through time vector
    for i=2:length( tvec )
        % Previous step information
        prev_state = integ_states(i-1,:)';
        prev_store = store_states(i-1,:)';
        
        [state_deriv, new_store] = feval( diff_eq, tvec(i-1), prev_state, prev_store );
        
        % Fill in info for this time step
        dt = tvec(i) - tvec(i-1);
        integ_states(i,:) = dt * state_deriv + prev_state;
        store_states(i,:) = new_store;
    end
    
    tinteg = tvec;
end