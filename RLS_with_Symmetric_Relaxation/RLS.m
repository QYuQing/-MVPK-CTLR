function [Y_pre, A] = RLS(Y, K, lambda)
    n = size(K, 1);
    A = (K + lambda * eye(n)) \ Y;
    Y_pre = K * A;
end