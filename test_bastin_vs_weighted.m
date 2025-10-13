function test_bastin_vs_weighted()
% Test script to verify that 'bastin' and 'weighted' methods give different results
% Uses same parameters as scan_SHC_MoverTz_DSM but with T=300K and fewer points for speed

%% ===== 物理常數與 SHC 設定 =====
a_angstrom = 1;
hbar = 1; e = 1;
Nk = 21;                   % Smaller grid for faster testing
eta_broad = 0.001;         % Kubo broadening (eV)
Ef   = 0.0; mu = 0.0; 
T = 300;                   % 300K to see difference between methods
kB = 8.617333262e-5;       % eV/K
T_eV = T * kB;             % Convert to eV
alpha = 'x';  beta = 'y';  gamma = 'z';   % σ^{s_z}_{xy}

%% ===== DSM 參數 =====
eta_vel = 0.89;
tz      = -3.4 * eta_vel;
beta4   = 0.67 * tz;
gamma4  = 0.335 * tz;

% Test with fewer points for speed
txy_over_tz_list = [1.0];              % Only one ratio
M_over_tz_grid   = linspace(0.5, 1.5, 10);  % Only 10 points

Sigma_bastin = zeros(1, numel(M_over_tz_grid));
Sigma_weighted = zeros(1, numel(M_over_tz_grid));

fprintf('Testing Bastin vs Weighted methods at T = %.0f K (%.4f eV)\n', T, T_eV);
fprintf('Using txy/tz = %.1f\n', txy_over_tz_list(1));
fprintf('M/tz range: [%.1f, %.1f] with %d points\n', ...
        M_over_tz_grid(1), M_over_tz_grid(end), numel(M_over_tz_grid));
fprintf('Progress: ');

%% ===== 主掃描 =====
txy = txy_over_tz_list(1) * tz;

for im = 1:numel(M_over_tz_grid)
    M_over_tz = M_over_tz_grid(im);
    Mval      = M_over_tz * tz;

    % --- 直接生成 ftn58sparse 兼容結構 ---
    ftn = build_ftn58sparse_DSM(eta_vel, txy, tz, Mval, beta4, gamma4);

    % --- SHC 前處理 ---
    params.ftn58 = ftn;
    params.Nk    = Nk;
    params.eta   = eta_broad;
    params.hbar  = hbar;
    params.electronic_charge = e;
    params.alpha = alpha; params.beta = beta; params.gamma = gamma;

    cache = shc.precompute_kgrid(params);

    % --- 評估 σ^{s_z}_{xy} using BASTIN method ---
    out_bastin = shc.eval_sigma(cache, mu, Ef, T_eV, 'bastin');
    Sigma_bastin(im) = real(out_bastin.sigma);
    
    % --- 評估 σ^{s_z}_{xy} using WEIGHTED method ---
    out_weighted = shc.eval_sigma(cache, mu, Ef, T_eV, 'weighted');
    Sigma_weighted(im) = real(out_weighted.sigma);
    
    fprintf('%d ', im);
    if mod(im, 10) == 0
        fprintf('\n         ');
    end
end
fprintf('\nDone!\n\n');

%% ===== 分析結果 =====
% Convert to plotting units
e2_over_h_S = -3.874045e-5;          % Siemens
a_meter = a_angstrom * 1e-10;       % m
scale = e2_over_h_S / a_meter;   % (Ω·m)^-1
Sigma_bastin_plot = scale * Sigma_bastin;
Sigma_weighted_plot = scale * Sigma_weighted;

% Calculate differences
abs_diff = abs(Sigma_bastin_plot - Sigma_weighted_plot);
rel_diff = abs_diff ./ max(abs(Sigma_bastin_plot), 1e-12) * 100; % Relative difference in %

fprintf('=== RESULTS ===\n');
fprintf('Max absolute difference: %.2e\n', max(abs_diff));
fprintf('Max relative difference: %.2f%%\n', max(rel_diff));
fprintf('Average relative difference: %.2f%%\n', mean(rel_diff));

if max(abs_diff) < 1e-10
    fprintf('WARNING: Methods give nearly identical results!\n');
    fprintf('This might indicate:\n');
    fprintf('  1. Implementation issue in eval_sigma\n');
    fprintf('  2. Physical regime where methods converge\n');
    fprintf('  3. Need higher temperature or different parameters\n');
else
    fprintf('SUCCESS: Methods give different results as expected!\n');
end

%% ===== 繪圖比較 =====
figure('Color','w', 'Position', [100, 100, 1200, 400]);

% Plot 1: Both methods
subplot(1,3,1);
plot(M_over_tz_grid, Sigma_bastin_plot, 'r-o', 'LineWidth', 2, 'DisplayName', 'Bastin');
hold on;
plot(M_over_tz_grid, Sigma_weighted_plot, 'b--s', 'LineWidth', 2, 'DisplayName', 'Weighted');
xlabel('M/t_z');
ylabel('\sigma^{\tilde z}_{xy} [(ħ/e)(\Omega\cdot m)^{-1}] \times 10^4');
title(sprintf('Both Methods (T=%.0fK)', T));
legend('Location', 'best');
grid on;

% Plot 2: Absolute difference
subplot(1,3,2);
semilogy(M_over_tz_grid, abs_diff, 'g-o', 'LineWidth', 2);
xlabel('M/t_z');
ylabel('|Bastin - Weighted|');
title('Absolute Difference');
grid on;

% Plot 3: Relative difference
subplot(1,3,3);
plot(M_over_tz_grid, rel_diff, 'm-o', 'LineWidth', 2);
xlabel('M/t_z');
ylabel('Relative Difference (%)');
title('Relative Difference');
grid on;

sgtitle(sprintf('Bastin vs Weighted Method Comparison (T=%.0fK, t_{xy}/t_z=%.1f)', ...
                T, txy_over_tz_list(1)));

%% ===== 保存結果 =====
outname = sprintf('test_bastin_vs_weighted_T%.0f', T);
savefig([outname '.fig']);
print(gcf, [outname '.jpg'], '-djpeg', '-r300');

fprintf('\nFigure saved as: %s.fig and %s.jpg\n', outname, outname);
end