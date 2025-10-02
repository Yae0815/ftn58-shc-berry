function sanity_check_DSM_full()
    % params
    eta=0.8; txy=0.6; tz=1.0; M=0.2; beta=0.3; gamma=0.1;
    ftn = build_ftn58sparse_DSM(eta,txy,tz,M,beta,gamma);

    % assemble H(k) from ftn (2π convention)
    function H = Hk(kx,ky,kz)
        H = zeros(4);
        phi = 2*pi*[kx,ky,kz];
        for n=1:size(ftn.tt,1)
            i = ftn.ij(n,1); j = ftn.ij(n,2);
            R = ftn.dd(n,:);
            H(i,j) = H(i,j) + ftn.tt(n)*exp(1i*dot(phi,R));
            if i~=j, H(j,i) = conj(H(i,j)); end
        end
    end

    % closed-form H(k) in basis [↑A, ↑B, ↓A, ↓B]
    BG = beta+gamma; BM = beta-gamma;
    function H = H_target(kx,ky,kz)
    cx = cos(2*pi*kx);  cy = cos(2*pi*ky);  cz = cos(2*pi*kz);
    sx = sin(2*pi*kx);  sy = sin(2*pi*ky);  sz = sin(2*pi*kz);

    BG = beta+gamma; BM = beta-gamma;

    % a3: σz ⊗ s0  （注意沒有 1/2！）
    a3 = M - txy*(cx+cy) - tz*cz;

    % a1: +eta sin kx (σx ⊗ sz)
    A1 = eta*sx;

    % a2: -eta sin ky (σy ⊗ s0)
    A2 = -eta*sy;

    % a4: (β+γ) sin kz (cos ky - cos kx) (σx ⊗ sx)
    A4 = BG*sz*(cy - cx);

    % a5: -(β-γ) sin kz sin kx sin ky (σx ⊗ sy)
    A5 = +BM*sz*sx*sy; 

    H = zeros(4);
    % a3
    H(1,1)=+a3; H(2,2)=-a3; H(3,3)=+a3; H(4,4)=-a3;

    % a1: σx ⊗ sz（↑:+σx, ↓:-σx）
    H(1,2)=H(1,2)+(+A1); H(2,1)=H(2,1)+(+A1);
    H(3,4)=H(3,4)+(-A1); H(4,3)=H(4,3)+(-A1);

    % a2: σy ⊗ s0（兩個自旋相同）
    H(1,2)=H(1,2)+(-1i*A2); H(2,1)=H(2,1)+(+1i*A2);
    H(3,4)=H(3,4)+(-1i*A2); H(4,3)=H(4,3)+(+1i*A2);

    % a4: σx ⊗ sx（↑A↔↓B, ↑B↔↓A）
    H(1,4)=H(1,4)+(+A4); H(4,1)=H(4,1)+(+A4);
    H(2,3)=H(2,3)+(+A4); H(3,2)=H(3,2)+(+A4);

    % a5: σx ⊗ sy（同樣配對，但帶 ±i）
    H(1,4)=H(1,4)+(+1i*A5); H(4,1)=H(4,1)+(-1i*A5);
    H(2,3)=H(2,3)+(+1i*A5); H(3,2)=H(3,2)+(-1i*A5);
end

    rng(0); errs = zeros(10,1); herm = zeros(10,1);
    for t=1:10
        k = rand(1,3);
        H1 = Hk(k(1),k(2),k(3));
        H2 = H_target(k(1),k(2),k(3));
        errs(t) = max(abs(H1(:)-H2(:)));
        herm(t) = norm(H1 - H1','fro');
    end
    fprintf('FULL max |H - H_target| = %.3e\n', max(errs));
    fprintf('Hermiticity check (max Fro norm) = %.3e\n', max(herm));
end
