%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Floquet multipliers are eigenvalues of a real mapping: this implies that
% they are either real or complex-conjugate pairs. In addition, α is ω-periodic.
% The two cases α = 0 and α = ω/2 are called harmonic and subharmonic,
% respectively, and correspond to positive or negative real Floquet multipliers,
% whereas intermediate values correspond to a complex Floquet multiplier.
%
% The procedure for the stability analysis is to fix the wavenumber 'k' and
% the exponents 'σ + iɑ', usually as σ=0 and at ɑ=0 or ɑ=1/2ω. We then solve 
% (3.15) for the eigenvalues 'a'. Only real and positive values of 'a' are 
% meaningful in this context.
%
function [results_harmonic, results_subharmonic] = sweep_over_omega(w_dim, N, k, nu, mu, rho, g, gamma, h, m_max)

  results_harmonic = [];
  results_subharmonic = [];

  for m = 1:m_max
    if (mod(m,100) == 0)
      m
    end

    omega = m * 60 * pi / m_max;
    
    sigma_values = 0;%linspace(0, 2, num_sigma_values) * g;

    % Harmonic case, α = 0
    alpha = 0;
    is_harmonic = true;
    a_harmonic = compute_eigenvalues(sigma_values, alpha, is_harmonic, k, N, omega, nu, mu, rho, g, gamma, h);

    % Sub-harmonic case, α = ω/2
    alpha = omega / 2;
    is_harmonic = false;
    a_subharmonic = compute_eigenvalues(sigma_values, alpha, is_harmonic, k, N, omega, nu, mu, rho, g, gamma, h);

    % Accumulate the results in a matrix
    results_harmonic = [results_harmonic; a_harmonic];
    results_subharmonic = [results_subharmonic; a_subharmonic];
  end
end