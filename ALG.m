%% Recursive algorithm for parameter estimation problem in linear regression model

    % Input: M - true covarience matrix
    % x0 - expection of vector x
    % H - observation matrix
    % Output: mr - root mean square error of estimation

function mr = ALG(x0, M, H)
    % Constants
    p = size(H, 1);         % iterations
    n = length(M);          % dimensions
    ns = 1000;              % number samples

    % Generate multivariate normal distribution
    x = mvnrnd(x0, M, ns)';

    % Calculate observation vector
    z = H*x;

    % Eigendecomposition or spectrum decomposition of M
    [U, D] = eig(M);

    % Resort eigenvalues descend and corresponding eigenvectors
    for i = 1:n-1
        for j = (i + 1):n
            if D(i, i) < D(j, j)
                tempvalue = D(i, i);
                D(i, i) = D(j, j);
                D(j, j) = tempvalue;
                tempvector = U(:, i);
                U(:, i) = U(:, j);
                U(:, j) = tempvector;
            end
        end
    end

    % Initialization
    for i = 1 : n
        M0{i} = D(1:i,1:i); U0{i} = U(:, 1:i);
        xe{i} = (U0{i}'*x0)*ones(1, ns);
    end
    
    M0{n+1} = eye(n); U0{n+1} = eye(n);
    xe{n+1} = (U0{n+1}'*x0)*ones(1, ns);

    % Main algorithm
    for i=1:p
        for j = 1: n+1
            h{j} = H(i,:)*U0{j};
            K{j} = M0{j}*h{j}'*MPInverse(h{j}*M0{j}*h{j}');
            M0{j} = M0{j} - K{j}*h{j}*M0{j};
            xe{j} = xe{j} + K{j}*(z(i, :) - h{j}*xe{j});
            xh{j} = U0{j}*xe{j};
        end
    end
    
    % Calculate error of estimation
    rms = zeros(n+1, ns);
    for i = 1: n
        for j = 1:n+1
           rms(j, :) =  rms(j, :) + (xh{j}(i,:) - x(i,:)).^2;
        end
    end
    rms = rms.^0.5;
    mr = mean(rms');
    
    % Plot root mean square error
    hold on;
    plot(rms(10,:), 'm');
    plot(rms(20,:), 'b');
    plot(rms(n,:), 'c--');
    plot(rms(n + 1,:), 'r:');
    legend('m = 10', 'm = 20', 'm = 100', 'Id');
    
    figure;
    i = p:n;
    plot(i, mr(p:n));
    ylabel('root mean square');
    xlabel('m');
end