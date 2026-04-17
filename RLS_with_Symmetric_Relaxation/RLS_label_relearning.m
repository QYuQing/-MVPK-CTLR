function [H, A, R, obj] = RLS_label_relearning(Y, K, lambda, maxIter, tol)
% RLS with label relearning
%
% Input:
%   Y       - n x C label matrix, usually one-hot
%   K       - n x n kernel matrix
%   lambda  - regularization parameter
%   maxIter - maximum number of alternating iterations
%   tol     - stopping tolerance
%
% Output:
%   H       - prediction matrix, H = K*A
%   A       - coefficient matrix
%   R       - relearned / enhanced label matrix
%   obj     - objective values

    if nargin < 4
        maxIter = 50;
    end
    if nargin < 5
        tol = 1e-6;
    end

    [n, C] = size(Y);

    % initialization
    R = Y;
    obj = zeros(maxIter, 1);

    I_n = eye(n);

    for t = 1:maxIter
        R_old = R;

        % ===== Step 1: optimize A with fixed R =====
        % min_A ||R - K*A||_F^2 + lambda * tr(A' K A)
        % closed-form:
        % A = (K + lambda*I)^(-1) R
        A = (K + lambda * I_n) \ R;
        H = K * A;

        % ===== Step 2: optimize R with fixed H =====
        for i = 1:n
            % true class index
            [~, yi] = max(Y(i, :));

            % update relearned label row r_i based on h_i
            R(i, :) = update_R_row(H(i, :), yi);
        end

        % ===== Objective =====
        obj(t) = norm(H - R, 'fro')^2 + lambda * trace(A' * K * A);

        % ===== Convergence check =====
        relchg = norm(R - R_old, 'fro') / max(1, norm(R_old, 'fro'));
        if relchg < tol
            obj = obj(1:t);
            fprintf('Converged at iteration %d, relchg = %.3e\n', t, relchg);
            break;
        end
    end

    if length(obj) == maxIter && obj(end) == 0
        obj = obj(1:t);
    end
end

function r = update_R_row(h, m)
% Update one row r_i from h_i according to Algorithm 2
%
% Input:
%   h - 1 x C prediction scores of sample i
%   m - true class index
%
% Output:
%   r - 1 x C relearned label row

    C = length(h);
    phi = -inf(1, C);

    % compute phi_j = h_ij + 1 - h_im for j ~= m
    for j = 1:C
        if j ~= m
            phi(j) = h(j) + 1 - h(m);
        end
    end

    zeta = 0;
    iter = 0;

    % Algorithm 2
    for j = 1:C
        if j ~= m
            if phi(j) > 0
                zeta = zeta + phi(j);
                iter = iter + 1;
            end
        end
    end

    zeta = zeta / (1 + iter);

    % update r with Eq. (34)
    r = zeros(1, C);
    for j = 1:C
        if j == m
            r(j) = h(j) + zeta;
        else
            r(j) = h(j) + min(zeta - phi(j), 0);
        end
    end
end