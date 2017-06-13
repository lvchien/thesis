%% Function generate a symmetric positive semi-definite matrix that be a covariance matrix

    % Input: n - demension of variables vector
    % Output: M - covariance matrix

function M = GenerateCovMatrix(n)
    A = randn(n);
    [U, ignore] = eig(A+A');
    M = U*diag(0.5*abs(randn(n,1)))*U';
end