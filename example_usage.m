% Example usage of the functions "nPSO_model" and "create_mixture_gaussian_gamma_pdf".
% The script generates and plots three examples:
% - nPSO model with 4 communities and default settings (Gaussian mixture distribution)
% - nPSO model with 8 communities and custom settings (Gaussian and Gamma mixture distribution)
% - PSO model

N = 100;
m = 4;
T = 0.1;
gamma = 3;
plot_flag = 1;

% nPSO model with 4 communities and default settings
distr = 4;
[x_nPSO1, coords_nPSO1, comm_nPSO1] = nPSO_model(N, m, T, gamma, distr, plot_flag);

% nPSO model with 8 communities and custom settings
distr = create_mixture_gaussian_gamma_pdf(8);
[x_nPSO2, coords_nPSO2, comm_nPSO2] = nPSO_model(N, m, T, gamma, distr, plot_flag);

% PSO model
distr = 0;
[x_PSO, coords_PSO] = nPSO_model(N, m, T, gamma, distr, plot_flag);