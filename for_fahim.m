% Main code to calculate eigenvalues for the parametric instability of the
% interface between two fluids using the approach to Kumar & Tuckerman (1994)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_dim = 512*pi*3.5e-3;                        % width, m
g_dim = 9.81;                                 % gravity, m/s²
rho_dim = [998.0, 1.25];                      % densities (bottom/top), kg/m³
mu_dim = [0.873e-3, 1.81e-5];                 % dyn. visc. Pa s
nu_dim = mu_dim./rho_dim;                     % kin. visc. m²/s
h_dim = [3.5e-3, 3.5e-3];                     % layer heights, m
gamma_dim = 0.0175;                           % interfacial tension, N/m²
omega_dim = 1142.86;                          % frequency, rad/s
f_dim = omega_dim/(2*pi);                     % frequency, Hz

m_max = 9600;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 20;                                       % Trunctation Fourier series
atwood = -diff(rho_dim)/sum(rho_dim);         % Atwood number

k0 = [1:m_max] * pi / w_dim;
omegan2 = atwood * g_dim .* k0 .* tanh(k0 * min(h_dim)) .* (1 - (gamma_dim * k0.^2) / diff(rho_dim));
omegan = sqrt(omegan2);               

[a_harmonic, a_subharmonic] = sweep_over_k(w_dim, N, omega_dim, nu_dim, mu_dim, rho_dim, g_dim, gamma_dim, h_dim, m_max);

figure(2)
k = a_harmonic(:,1);
plot(k, a_harmonic(:,3), '.k')
hold on 
plot(k, a_subharmonic(:,3), '.r')
for i = 4:12
  plot(k, a_harmonic(:,i), '.k')
  plot(k, a_subharmonic(:,i), '.r')
end
hold off
xlabel('$k$ (m$^{-1}$)', 'Interpreter', 'latex')
ylabel('$a/g$', 'Interpreter', 'latex')

% Add the vertical line at x = 2*pi/w_dim
xline(2*pi/3.5e-3, '--b', 'Label', '$k=2\pi/w$','Interpreter','latex');


% Set figure properties
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 4 3]); % Set the figure size (width x height)
set(gcf, 'PaperSize', [4 3]); % Ensure the paper size matches the position

% Set font size and line width
set(gcf, 'DefaultAxesFontSize', 12);
set(gcf, 'DefaultTextFontSize', 12);
set(gca, 'LineWidth', 1.5);

% Set the resolution (DPI)
set(gcf, 'PaperPositionMode', 'auto');
print('-depsc2', '-r300', 'for_fahim.eps'); % Save as EPS with 300 DPI
close all

save('for_fahim.mat', "a_harmonic", "a_subharmonic")