%%%
% File: SpineLabParticleStep.m
% Author: Calvin Kuo
% Date: 11-03-2018

function [x_t_t, sig_t_t, new_particle_states] = SpineLabParticleStep( x_t1_t1, sig_t1_t1, particle_states, alpha, params, dt, tnow )

    % Extract parameters
    tc = params.tc;
    tn_list = params.tn_list;
    
    sigAlpha = params.sigAlpha;
    sigAfferent_list = params.sigAfferent_list; % For separating out affernent noises, regular vs. irregular
    sigPrior = params.sigPrior;
    
    Tcan = params.Tcan;
    
    % Build Particle Update Matrix A (state dynamics), B (input transform),
    % C (noise transform)
    nStates = size( particle_states, 1 );
    assert( mod( nStates, 6 ) == 0 );
    nLines = nStates / 6;
    A = zeros( nStates, nStates );
    B = zeros( nStates, 3 );
    C = zeros( nStates, nStates );
    F = zeros( nStates, nStates );
    FnoiseVec = zeros( nStates, 1 ); 
    noiseVec = zeros( nStates, 1 );
    for i=1:nLines
        Coff = (i-1)*3;
        Doff = (i-1)*3 + nLines*3;
        A(Coff+1:Coff+3, Coff+1:Coff+3) = eye(3)*(1-dt/tc);
        A(Doff+1:Doff+3, Coff+1:Coff+3) = eye(3)*(-dt/tc);
        A(Doff+1:Doff+3, Doff+1:Doff+3) = eye(3)*(1-dt/tn_list(i));
        
        B(Coff+1:Coff+3,:) = eye(3)*dt;
        B(Doff+1:Doff+3,:) = eye(3)*dt;
        
        C(Coff+1:Coff+3, Coff+1:Coff+3) = eye(3)*dt;
        C(Doff+1:Doff+3, Coff+1:Coff+3) = eye(3)*dt;
        F(Doff+1:Doff+3, Doff+1:Doff+3) = eye(3);
        
        noiseVec(Coff+1:Coff+3) = sqrt( sigAlpha ) * randn(3,1);
        FnoiseVec(Doff+1:Doff+3) = sqrt( sigAfferent_list(:,:,i) ) * randn(3,1);
    end
    
    % CANAL AND LONG TERM ADAPTATION STATE DYNAMICS
    new_particle_states = A * particle_states + B * Tcan * alpha + C * noiseVec;
    afferent_states = new_particle_states + F * FnoiseVec;
    %afferent_states = new_particle_states;
    % Particle weighting
    kk = reshape( afferent_states((nStates/2+1):end), 3, nLines );
    compare_means = [ 0,0,0; x_t1_t1' ];% [x_t1_t1'];%[ 0,0,0; x_t1_t1' ];
    compare_sigs = [diag(sigPrior)'; diag( sig_t1_t1 )']; %, [0.28, 0.28, 0.28] )]; %[diag(sigPrior)'; diag( sig_t1_t1 )' + [0.01, 0.01, 0.01]]; %[diag(sigPrior)'; diag( sig_t1_t1 )' + [0.01, 0.01, 0.01]]; %;
    [x_t_t, sig_t_t_vec, wts] = CanalParticleWeighting( kk, compare_means, compare_sigs );
    sig_t_t = diag( sig_t_t_vec );
    
    ll = Tcan'*kk;
    x_rot = Tcan'*x_t_t;
    sig_rot = Tcan'*sig_t_t_vec;
    plot( ones(nLines,1)*tnow, ll(2,:), 'x', 'Color', [0.7, 0.7, 0.7],'HandleVisibility','off' )
    plot( tnow, mean( ll(2,:) ), 'ro','HandleVisibility','off' );
    plot( tnow, x_rot(2), 'gs', 'MarkerSize', 10.0,'HandleVisibility','off' );
    plot( [tnow, tnow], [x_rot(2) - sig_rot(2), x_rot(2) + sig_rot(2)], 'k-', 'LineWidth', 2,'HandleVisibility','off'  )
end