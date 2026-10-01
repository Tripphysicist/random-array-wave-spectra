function DS = bdm_nnls_method(Sxyn, Gwt, theta, fi, k, opt)
% BDM_NNLS_METHOD Bayesian Directional Spectrum via Non-Negative Least Squares
%
% This is the rigorous port of the NOWPHAS ABIC-NNLS operational algorithm.
% It uses True ABIC (with log-determinant) and a two-pass hyperparameter search.
[m, nt, nf] = size(Gwt);
DS = repmat(1/(2*pi), nt, nf);
dtheta = theta(2) - theta(1);
% 1. Create the Circular 2nd-Derivative Smoothing Matrix (D_mat)
D_mat = diag(-2*ones(nt,1)) + diag(ones(nt-1,1), 1) + diag(ones(nt-1,1), -1);
D_mat(1, nt) = 1;
D_mat(nt, 1) = 1;
for ff = k
    % --- Extract Valid Equations ---
    M = 0;
    % FIX: Initialize as correctly dimensioned empty arrays
    A = zeros(0, nt); 
    B = zeros(0, 1);
    tol = sqrt(eps);
    
    for ix = 1:m
        for iy = ix:m
            Htemp = Gwt(ix,:,ff) .* conj(Gwt(iy,:,ff));
            if any(abs(diff(real(Htemp))) > tol)
                M = M + 1;
                A(M, :) = real(Htemp) * dtheta;
                B(M, 1) = real(Sxyn(ix, iy, ff));
            end
            if any(abs(diff(imag(Htemp))) > tol)
                M = M + 1;
                A(M, :) = imag(Htemp) * dtheta;
                B(M, 1) = imag(Sxyn(ix, iy, ff));
            end
        end
    end
    
    % FIX: If no valid equations exist (e.g. at f=0), skip the solver.
    % The spectrum will default to the uniform distribution initialized above.
    if M == 0
        continue;
    end
    
    % --- Normalization ---
    % Scale B to prevent floating point underflow during residual squaring
    scale_factor = max(abs(B));
    if scale_factor > 0
        B = B / scale_factor;
    end
    
    % --- Two-Pass Hyperparameter Search ---
    % Pass 1: Coarse Search
    u_list_coarse = logspace(-4, 2, 20);
    [best_u, ~] = run_abic_sweep(A, B, D_mat, u_list_coarse, nt, M);
    
    % Pass 2: Fine Search (Mimics Golden Section)
    u_list_fine = linspace(best_u * 0.1, best_u * 10, 20);
    [~, best_x] = run_abic_sweep(A, B, D_mat, u_list_fine, nt, M);
    
    % De-normalize and store
    DS(:, ff) = best_x * scale_factor;
end
% Final WAFO normalization
DS = normspfn(DS, theta);
end
% --- Helper Function for ABIC Optimization ---
function [best_u, best_x] = run_abic_sweep(A, B, D_mat, u_list, K, N_eq)
    best_ABIC = inf;
    best_u = u_list(1);
    best_x = zeros(K, 1);
    
    % Silence the lsqnonneg solver
    opts = optimset('Display', 'none');
    AtA = A' * A;
    DtD = D_mat' * D_mat;
    
    for u = u_list
        u2 = u^2;
        % Stack matrices for NNLS
        C = [A; u * D_mat];
        d = [B; zeros(K, 1)];
        
        % Solve bounded linear system silently
        [x, ~, ~, ~, ~] = lsqnonneg(C, d, opts);
        
        % --- True ABIC Calculation ---
        % 1. Variance of the residuals
        sq_sigma = (norm(A*x - B)^2 + u2 * norm(D_mat*x)^2) / N_eq;
        
        if sq_sigma <= 0
            continue;
        end
        
        % 2. Log-Determinant of the Hessian (Symmetric Positive Semi-Definite)
        Hessian = AtA + u2 * DtD;
        s_vals = eig((Hessian + Hessian') / 2);
        ln_det = sum(log(s_vals(s_vals > eps))); 
        
        % 3. Final ABIC score
        ABIC = N_eq * log(sq_sigma) - K * log(u2) + ln_det;
        
        if ABIC < best_ABIC
            best_ABIC = ABIC;
            best_u = u;
            best_x = x;
        end
    end
end