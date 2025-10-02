function scan_SHC_MoverTz_DSM()
% 掃描 DSM：固定 tz，txy/tz ∈ {0.5,1.0,1.5}，M/tz ∈ [0,2]
% 直接用 build_ftn58sparse_DSM() 回傳的 struct 進 +shc，不再跑 TBHmftn。

clc; clear;

%% ===== 物理常數與 SHC 設定 =====
hbar = 6.582119569e-16;   % eV·s
e    = 1.602176634e-19;   % C
Nk   = 21;                % k-mesh
eta_broad = 0.001;         % Kubo broadening (eV)
Ef   = 0.0; mu = 0.0; T = 0;      % Fermi level/ref
alpha = 'x';  beta = 'y';  gamma = 'z';   % σ^{s_z}_{xy}

%% ===== DSM 參數 =====
eta_vel = 0.89;
tz      = -3.4 * eta_vel;
beta4   = 0.67 * tz;
gamma4  = 0.335 * tz;

txy_over_tz_list = [0.5, 1.0, 1.5];
M_over_tz_grid   = linspace(0, 2.0, 11);

Sigma = zeros(numel(txy_over_tz_list), numel(M_over_tz_grid));

%% ===== 主掃描 =====
for ir = 1:numel(txy_over_tz_list)
    ratio = txy_over_tz_list(ir);
    txy   = ratio * tz;

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

        % --- 評估 σ^{s_z}_{xy}（Kubo-Berry weighted）---
        out = shc.eval_sigma(cache, mu, Ef, T, 'weighted');
        Sigma(ir, im) = real(out.sigma);
    end
end

%% ===== 繪圖 =====
Sigma_plot = -1e4 * (hbar/e) * Sigma;    % 轉成 [(ħ/e)(Ω·m)^{-1}] ×10^4

figure('Color','w'); hold on;
for iRatio = 1:numel(txy_over_tz_list)
    plot(M_over_tz_grid, Sigma_plot(iRatio,:), 'LineWidth', 2);
end
xlabel('M / t_z','FontSize',12,'Interpreter','latex');

ylabel('$\sigma^{\tilde z}_{xy}\;[(\hbar/e)(\Omega\!\cdot\! m)^{-1}]\times 10^{4}$',...
       'FontSize',12,'Interpreter','latex');

legend(arrayfun(@(r) sprintf('t_{xy}/t_z = %.1f', r), ...
       txy_over_tz_list, 'UniformOutput', false), 'Location','best');

grid on; box on;
title(sprintf('DSM: \\eta=%.2f, t_z=%.3f eV, Nk=%d, \\eta_{broad}=%.3f eV', ...
      eta_vel, tz, Nk, eta_broad), 'Interpreter','tex');