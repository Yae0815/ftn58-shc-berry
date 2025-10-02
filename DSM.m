%% ===== DSM ftn58 generator (spin ⊗ orbit; 2π phase) =====
% Basis (spin ⊗ orbit): 1=|↑A>, 2=|↑B>, 3=|↓A>, 4=|↓B>
% Convention: H(k) = sum_R H(R) * exp(i * 2π * k · R)

clear; clc;

%% ===== Parameters for DSM =====
eta   = 0.89;           % velocity scale (eV)
tz    = -3.4 * eta;     % out-of-plane mass (eV)
txy   = 0.3;            % in-plane mass (eV)
M     = -0.1;           % Dirac mass (eV)
beta  = 0.67  * tz;     % for a4, a5 (eV)
gamma = 0.335 * tz;     % for a4, a5 (eV)
BG    = beta + gamma;
BM    = beta - gamma;

%% ===== a3: σz ⊗ s0 (onsite + cos) =====
% Note: include BOTH +R and -R so diagonal stays real (cos) and H is Hermitian.

ftn58_a3_diag = [
    1 1  (+M)   0 0 0
    2 2  (-M)   0 0 0
    3 3  (+M)   0 0 0
    4 4  (-M)   0 0 0
];

ftn58_a3_cos = [
    % ±x : -txy * cos kx
    1 1 (-txy/2)  +1 0 0
    2 2 (+txy/2)  +1 0 0
    3 3 (-txy/2)  +1 0 0
    4 4 (+txy/2)  +1 0 0
    1 1 (-txy/2)  -1 0 0
    2 2 (+txy/2)  -1 0 0
    3 3 (-txy/2)  -1 0 0
    4 4 (+txy/2)  -1 0 0

    % ±y : -txy * cos ky
    1 1 (-txy/2)   0 +1 0
    2 2 (+txy/2)   0 +1 0
    3 3 (-txy/2)   0 +1 0
    4 4 (+txy/2)   0 +1 0
    1 1 (-txy/2)   0 -1 0
    2 2 (+txy/2)   0 -1 0
    3 3 (-txy/2)   0 -1 0
    4 4 (+txy/2)   0 -1 0

    % ±z : -tz  * cos kz
    1 1 (-tz/2)    0 0 +1
    2 2 (+tz/2)    0 0 +1
    3 3 (-tz/2)    0 0 +1
    4 4 (+tz/2)    0 0 +1
    1 1 (-tz/2)    0 0 -1
    2 2 (+tz/2)    0 0 -1
    3 3 (-tz/2)    0 0 -1
    4 4 (+tz/2)    0 0 -1
];

%% ===== a1: +eta * sin kx * (σx ⊗ s_z) =====
% Put ±x with opposite signs so it sums to sin(2π kx).
ftn58_a1 = [
    % +x
    1 2 ( +eta/(2i))  +1 0 0
    3 4 ( -eta/(2i))  +1 0 0
    % -x
    1 2 ( -eta/(2i))  -1 0 0
    3 4 ( +eta/(2i))  -1 0 0
];

%% ===== a2: -eta * sin ky * (σy ⊗ s0) =====
% Put ±y with opposite signs so it sums to sin(2π ky) (pure ±i off-diagonal).
ftn58_a2 = [
    % +y
    1 2 ( +eta/2)   0 +1 0
    3 4 ( +eta/2)   0 +1 0
    % -y
    1 2 ( -eta/2)   0 -1 0
    3 4 ( -eta/2)   0 -1 0
];

%% ===== a4: (β+γ) sin kz (cos ky - cos kx) * (s_x ⊗ σ_x) =====
% Use pure-imag coeffs and include z=±1 so only sin kz survives.
ftn58_a4 = [
    % +z, +ky part:  + BG/(4i) * e^{i2π(kz+ky)}
    1 4 ( +BG/(4i))   0 +1 +1
    2 3 ( +BG/(4i))   0 +1 +1
    % +z, -ky part
    1 4 ( +BG/(4i))   0 -1 +1
    2 3 ( +BG/(4i))   0 -1 +1
    % +z, -kx part:  - BG/(4i) * e^{i2π(kz+kx)}
    1 4 ( -BG/(4i))  +1  0 +1
    2 3 ( -BG/(4i))  +1  0 +1
    1 4 ( -BG/(4i))  -1  0 +1
    2 3 ( -BG/(4i))  -1  0 +1

    % -z partners with opposite overall signs to produce sin kz
    % ky block (flip sign)
    1 4 ( -BG/(4i))   0 +1 -1
    2 3 ( -BG/(4i))   0 +1 -1
    1 4 ( -BG/(4i))   0 -1 -1
    2 3 ( -BG/(4i))   0 -1 -1
    % kx block (flip sign of x-block relative to +z entries)
    1 4 ( +BG/(4i))  +1  0 -1
    2 3 ( +BG/(4i))  +1  0 -1
    1 4 ( +BG/(4i))  -1  0 -1
    2 3 ( +BG/(4i))  -1  0 -1
];

%% ===== a5: +(β−γ) sin kz sin kx sin ky * (s_y ⊗ σ_x) =====
% Build with real coeffs at z=±1 and opposite signs so it becomes i * sin kz * sin kx * sin ky.
ftn58_a5 = [
    % z=+1
    1 4 ( -BM/8)  +1 +1 +1
    2 3 ( -BM/8)  +1 +1 +1
    1 4 ( +BM/8)  -1 +1 +1
    2 3 ( +BM/8)  -1 +1 +1
    1 4 ( +BM/8)  +1 -1 +1
    2 3 ( +BM/8)  +1 -1 +1
    1 4 ( -BM/8)  -1 -1 +1
    2 3 ( -BM/8)  -1 -1 +1
    % z=-1 (opposite signs)
    1 4 ( +BM/8)  +1 +1 -1
    2 3 ( +BM/8)  +1 +1 -1
    1 4 ( -BM/8)  -1 +1 -1
    2 3 ( -BM/8)  -1 +1 -1
    1 4 ( -BM/8)  +1 -1 -1
    2 3 ( -BM/8)  +1 -1 -1
    1 4 ( +BM/8)  -1 -1 -1
    2 3 ( +BM/8)  -1 -1 -1
];

%% ----- collect & write ftn58.mat -----
ftn58_min = [ftn58_a3_diag; ftn58_a3_cos; ftn58_a1; ftn58_a2];
ftn58     = [ftn58_min; ftn58_a4; ftn58_a5];

nbond = size(ftn58,1);
ftn58  = [(1:nbond)' ftn58];                % prepend bond index
ftn58  = [[4 nbond 0 0 0 0 0]; ftn58];      % header: [norb, nbond, ... zeros]
save ftn58.mat ftn58;

% (optional) quick run
% eval('run TBHmftn.m');
% eval('run band_ftn.m');
