%% ===== Parameters for DSM (use these; old A1..D2 not used below) =====
eta = 0.89;    % velocity scale (eV)
tz  = -3.4 * eta;    % out-of-plane mass term (eV)
txy = 0.3;    % in-plane mass term (eV)
M   = -0.1;   % Dirac mass (eV)
beta = 0.67 * tz;  % for a4, a5 term (eV)
gamma= 0.335 * tz;  % for a4 ,a5 term (eV)

%% ===== a3: sigma_z s0 (onsite + cos terms) =====
% Onsite: M * sigma_z s0  => diag(+M,-M,+M,-M) at R=(0,0,0)
ftn58_a3_diag = [
    1 1  (+M)   0 0 0
    2 2  (-M)   0 0 0
    3 3  (+M)   0 0 0
    4 4  (-M)   0 0 0
];

% -txy * cos kx -> R=+x with coeff (-txy/2) * sigma_z s0 (diagonal)
% -txy * cos ky -> R=+y with coeff (-txy/2) * sigma_z s0 (diagonal)
% -tz  * cos kz -> R=+z with coeff (-tz /2) * sigma_z s0 (diagonal)
ftn58_a3_cos = [
    % R = +x
    1 1  (-txy/2)   1 0 0
    2 2  (+txy/2)   1 0 0
    3 3  (-txy/2)   1 0 0
    4 4  (+txy/2)   1 0 0
    % R = +y
    1 1  (-txy/2)   0 1 0
    2 2  (+txy/2)   0 1 0
    3 3  (-txy/2)   0 1 0
    4 4  (+txy/2)   0 1 0
    % R = +z
    1 1  (-tz/2)    0 0 1
    2 2  (+tz/2)    0 0 1
    3 3  (-tz/2)    0 0 1
    4 4  (+tz/2)    0 0 1
];

%% ===== a1: +eta * sin kx * (sigma_x s_z) =====
% Only upper-right (i<j) and R=+x; Hermitian/(-R)由腳本自動補
% For spin up (1<->2): +eta/(2i); for spin down (3<->4): -eta/(2i)
ftn58_a1 = [
    1 2  ( eta/(2i))   1 0 0   % up block
    3 4  (-eta/(2i))   1 0 0   % down block
];

%% ===== a2: -eta * sin ky * (sigma_y s0) =====
% After combining sigma_y and sin ky, the hopping is REAL:
% +y bond gives +eta/2 (and -y would be -eta/2, but we only list +y)
ftn58_a2 = [
    1 2  (+eta/2)   0 1 0
    3 4  (+eta/2)   0 1 0
];

%% ===== a4 term: (beta+gamma) * sin(kz) * (cos(ky) - cos(kx)) * (s_x ⊗ sigma_x)
% Define parameter for convenience:
BG = beta + gamma;   % set this to (beta + gamma) in eV, e.g., BG = beta + gamma;

% Coeff reminder:  +sin kz cos ky ⇒  +BG/(4i) at R=(0, ±1, +1)
%                   -sin kz cos kx ⇒  -BG/(4i) at R=(±1, 0, +1)
% Only upper-right entries (i<j): channels (1,4) and (2,3)

ftn58_a4 = [
    % --- from  sin(kz) cos(ky): R = (0, +1, +1) and (0, -1, +1)
    1 4  ( +BG/(4i) )   0  1  1
    2 3  ( +BG/(4i) )   0  1  1
    1 4  ( +BG/(4i) )   0 -1  1
    2 3  ( +BG/(4i) )   0 -1  1

    % --- from -sin(kz) cos(kx): R = (+1, 0, +1) and (-1, 0, +1)
    1 4  ( -BG/(4i) )   1  0  1
    2 3  ( -BG/(4i) )   1  0  1
    1 4  ( -BG/(4i) )  -1  0  1
    2 3  ( -BG/(4i) )  -1  0  1
];

%% ===== a5 term: -(beta - gamma) * sin(kz) * sin(kx) * sin(ky) * (s_y ⊗ sigma_x)
% Only list z=+1; your script will auto-hermit to add the -z partners.
% Channels (upper-right only): (1,4) and (2,3)
% Coefficient rule for z=+1:  t(R) = -(BM)/8 * (sx * sy), where R=(sx, sy, +1)

BM = beta - gamma;

ftn58_a5 = [
    % R = (+1, +1, +1):  sx*sy = +1  ->  t = -(BM)/8
    1 4  ( -BM/8 )    1  1  1
    2 3  ( -BM/8 )    1  1  1

    % R = (-1, +1, +1):  sx*sy = -1  ->  t = +(BM)/8
    1 4  ( +BM/8 )   -1  1  1
    2 3  ( +BM/8 )   -1  1  1

    % R = (+1, -1, +1):  sx*sy = -1  ->  t = +(BM)/8
    1 4  ( +BM/8 )    1 -1  1
    2 3  ( +BM/8 )    1 -1  1

    % R = (-1, -1, +1):  sx*sy = +1  ->  t = -(BM)/8
    1 4  ( -BM/8 )   -1 -1  1
    2 3  ( -BM/8 )   -1 -1  1
];


%% ----- collect minimal DSM terms -----
ftn58_min = [ftn58_a3_diag; ftn58_a3_cos; ftn58_a1; ftn58_a2];

% ----------------------------------------------------------------------------------
ftn58 = [ftn58_min; ftn58_a4; ftn58_a5];                 % 用最小DSM a1~a3
nbond = size(ftn58,1);
ftn58  = [(1:nbond)' ftn58];
ftn58  = [[4 nbond 0 0 0 0 0]; ftn58];

save ftn58.mat ftn58
eval('run TBHmftn.m');
eval('run band_ftn.m');
