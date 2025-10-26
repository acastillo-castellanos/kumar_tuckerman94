%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The coefficients A, and hence the system (3.12) depend on the Floquet exponent
% σ + iɑ via (3.8) and the boundary conditions (see Appendix B) in a complicated
% manner. However, the amplitude a of the external forcing appears linearly. 
% In fact, for fixed Floquet exponent, (3.12) can be considered to be an 
% eigenvalue problem with eigenvalues 'a' and eigenvectors 'v' whose components 
% are the real and imaginary parts of the ζ_n. 
% That is, we write (3.12) as the generalized eigenvalue problem
%
%   A ζ_n = a B ζ_n
%
% In (3.15), A is a diagonal complex matrix (see matrix_An) and B is a banded 
% matrix whose structure depends on ɑ (see banded_matrix_B).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a] = compute_eigenvalues(sigma_values, alpha, is_harmonic, k, N, omega, nu, mu, rho, g, gamma, h)
  num_sigma_values = length(sigma_values);
  An_block = cell(N + 1, 1);
  a = nan(num_sigma_values, 13);

  % Sweep over temporal frequencies
  for ii = 1:num_sigma_values
    sigma = sigma_values(ii);
    [B] = banded_matrix_B(N, is_harmonic);

    % Sweep over Fourier modes frequencies
    for ij = 1:N + 1
      n = ij - 1;
      An = matrix_An(sigma, alpha, n, omega, k, nu, mu, rho, g, gamma, h);
      An_block{ij} = [real(An) -imag(An); imag(An) real(An)];
    end

    % Assemble the block diagonal matrix A using the matrices in the cell array D
    A = blkdiag(An_block{:});

    % Solve the generalized eigenvalue problem A * x = lambda * (1/2*(rho(2)-rho(1))*(k^2)*B) * x
    [~, evals] = eig(A, 1 / 2 * (rho(2) - rho(1)) * (k^2) * g * full(B), "vector");
    

    % Filter eigenvalues to keep only those with:
    % - Imaginary parts close to zero
    % - Real parts that are positive and less than 2*g
    evals = evals(isinf(evals) == false);
    evals = evals(abs(imag(evals)) < 1e-12);
    evals = evals(abs(real(evals)) < 5.0);

    % Sort the remaining eigenvalues in ascending order
    evals = sort(evals(evals > 0));

    nevals = length(evals);
    a(ii,1) = k;
    a(ii,2) = omega;
    a(ii,3) = sigma;
    a(ii,4:3+min(nevals,10)) = evals(1:min(nevals,10));
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%