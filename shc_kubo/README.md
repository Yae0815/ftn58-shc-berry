# shc_kubo (Spin Hall via Kubo, μ-diagram)

## Install
1. 將整個 `shc_kubo/` 夾到任意位置。
2. 在 MATLAB 內 `addpath(genpath('<shc_kubo 路徑>'))` 或跑 `demo/demo_mu_diagram.m`。

## Inputs
- 需要一個 `ftn58sparse.mat`（含 fields: `ij, tt, dd, norb`），放在 `demo/` 的上一層（或自行修改路徑）。
- 所有 k 皆以 reduced coords 表示（phase = exp(i 2π k·R)）。

## Usage
- 修改 `demo/demo_mu_diagram.m` 的 `params`（αβγ、Nk、η、T、method）。
- 執行後自動掃 μ（`mu_grid`）並繪圖。

## Notes
- 單位為「晶格單位」下的 \(\sigma\)，方便觀察趨勢與符號；要絕對 SI 值需再處理體積與 \((2\pi)^3\) 因子。
- 避免 nested parfor：外層 μ 掃描用 for；內層 `spin_hall_main` 對 kx 用 parfor。
- 若 Γ 點 `[S^γ, H]` Frobenius norm 很大，表示自旋非守恆，請固定同一套 spin current 定義做比較。
