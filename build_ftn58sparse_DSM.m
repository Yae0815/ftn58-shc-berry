function ftn = build_ftn58sparse_DSM(eta, txy, tz, M, beta, gamma)
% Return ftn struct compatible with +shc/make_builders:
%   ftn.norb = 4
%   ftn.ij   = [nbond x 2] (i,j)   with only upper-tri + diagonals
%   ftn.tt   = [nbond x 1] complex hopping amplitudes
%   ftn.dd   = [nbond x 3] integer R-vectors (lattice units)
%
% Basis (spin ⊗ orbit, spin在前):
%   1=|↑A>, 2=|↑B>, 3=|↓A>, 4=|↓B>
% Convention: H(k) = sum_R H(R) e^{ i 2π k·R }  (your +shc uses 2π)

% ----- a3: diag + cos -----
a3_diag = [
    1 1  (+M)   0 0 0
    2 2  (-M)   0 0 0
    3 3  (+M)   0 0 0
    4 4  (-M)   0 0 0
];
a3_cos = [
    % +x
    1 1  (-txy/2)   1  0  0
    2 2  (+txy/2)   1  0  0
    3 3  (-txy/2)   1  0  0
    4 4  (+txy/2)   1  0  0
    % -x  (add these)
    1 1  (-txy/2)  -1  0  0
    2 2  (+txy/2)  -1  0  0
    3 3  (-txy/2)  -1  0  0
    4 4  (+txy/2)  -1  0  0

    % +y
    1 1  (-txy/2)   0  1  0
    2 2  (+txy/2)   0  1  0
    3 3  (-txy/2)   0  1  0
    4 4  (+txy/2)   0  1  0
    % -y  (add these)
    1 1  (-txy/2)   0 -1  0
    2 2  (+txy/2)   0 -1  0
    3 3  (-txy/2)   0 -1  0
    4 4  (+txy/2)   0 -1  0

    % +z
    1 1  (-tz/2)    0  0  1
    2 2  (+tz/2)    0  0  1
    3 3  (-tz/2)    0  0  1
    4 4  (+tz/2)    0  0  1
    % -z  (add these)
    1 1  (-tz/2)    0  0 -1
    2 2  (+tz/2)    0  0 -1
    3 3  (-tz/2)    0  0 -1
    4 4  (+tz/2)    0  0 -1
];

% ----- a1: +eta sin kx * (σx s_z) -----
a1 = [
    % +x
    1 2  ( +eta/(2i))   1  0  0
    3 4  ( -eta/(2i))   1  0  0
    % -x  (opposite sign to kill cos and keep sin)
    1 2  ( -eta/(2i))  -1  0  0
    3 4  ( +eta/(2i))  -1  0  0
];

% ----- a2: -eta sin ky * (σy s0) -----
a2 = [
    % +y
    1 2  (+eta/2)   0  1  0
    3 4  (+eta/2)   0  1  0
    % -y (opposite sign to kill cos and keep sin)
    1 2  (-eta/2)   0 -1  0
    3 4  (-eta/2)   0 -1  0
];


% ----- a4: (beta+gamma) sin kz (cos ky - cos kx) * (s_x ⊗ σ_x) -----
BG = beta + gamma;
a4 = [
    1 4  ( +BG/(4i) )   0  1  1
    2 3  ( +BG/(4i) )   0  1  1
    1 4  ( +BG/(4i) )   0 -1  1
    2 3  ( +BG/(4i) )   0 -1  1
    1 4  ( -BG/(4i) )   1  0  1
    2 3  ( -BG/(4i) )   1  0  1
    1 4  ( -BG/(4i) )  -1  0  1
    2 3  ( -BG/(4i) )  -1  0  1
    % ky terms at z = -1 (opposite sign)
    1 4  ( -BG/(4i) )   0  1 -1
    2 3  ( -BG/(4i) )   0  1 -1
    1 4  ( -BG/(4i) )   0 -1 -1
    2 3  ( -BG/(4i) )   0 -1 -1
    % kx terms at z = -1 (opposite sign of the x-block, so flip the sign)
    1 4  ( +BG/(4i) )   1  0 -1
    2 3  ( +BG/(4i) )   1  0 -1
    1 4  ( +BG/(4i) )  -1  0 -1
    2 3  ( +BG/(4i) )  -1  0 -1
];

% ----- a5: -(beta-gamma) sin kz sin kx sin ky * (s_y ⊗ σ_x) -----
BM = beta - gamma;
a5 = [
    1 4  ( -BM/8 )    1  1  1
    2 3  ( -BM/8 )    1  1  1
    1 4  ( +BM/8 )   -1  1  1
    2 3  ( +BM/8 )   -1  1  1
    1 4  ( +BM/8 )    1 -1  1
    2 3  ( +BM/8 )    1 -1  1
    1 4  ( -BM/8 )   -1 -1  1
    2 3  ( -BM/8 )   -1 -1  1
    1 4  ( +BM/8 )    1  1 -1
    2 3  ( +BM/8 )    1  1 -1
    1 4  ( -BM/8 )   -1  1 -1
    2 3  ( -BM/8 )   -1  1 -1
    1 4  ( -BM/8 )    1 -1 -1
    2 3  ( -BM/8 )    1 -1 -1
    1 4  ( +BM/8 )   -1 -1 -1
    2 3  ( +BM/8 )   -1 -1 -1
];

core = [a3_diag; a3_cos; a1; a2; a4; a5];  % [i j val a b c]

ftn.norb = 4;
ftn.ij   = core(:,1:2);
ftn.tt   = core(:,3);
ftn.dd   = core(:,4:6);

end
