%% MX-ONMTF Demo: Real Multiplex Network
% This script demonstrates how to run MX-ONMTF on real multiplex networks
% where ground truth community labels are not available.
%
% Reference:
%   M. Ortiz-Bouza and S. Aviyente, "Community Detection in Multiplex
%   Networks Based on Orthogonal Nonnegative Matrix Tri-Factorization,"
%   IEEE Access, vol. 12, pp. 6423-6436, 2024.

%% Add paths
addpath('..');
addpath('../helpers');

%% Step 0: Load your data
% Replace this section with your own data loading.
% Al should be a 1xL cell array where Al{l} is the n x n adjacency matrix
% of layer l.
%
% Example:
%   Al = {A1, A2, A3};  % 3-layer multiplex network
%
% Requirements:
%   - Each Al{l} must be square and symmetric
%   - All layers must have the same number of nodes n

% --- Placeholder: generate example data ---
n = 60;
L = 2;
GT1 = [ones(20,1); 2*ones(20,1); 3*ones(20,1)];
GT2 = [ones(20,1); 2*ones(20,1); 4*ones(20,1)];
A1 = generate_planted_partition(GT1, 0.4, 0.05);
A2 = generate_planted_partition(GT2, 0.4, 0.05);
Al = {A1, A2};
% --- End placeholder ---

%% Step 1: Find the number of communities
% This step is required for real data since we don't know the community
% structure in advance.

[kc, k, kpl] = findingk(Al);

fprintf('Detected communities:\n');
fprintf('  Common communities (kc): %d\n', kc);
fprintf('  Total per layer (k):     [%s]\n', num2str(k));
fprintf('  Private per layer (kpl): [%s]\n', num2str(kpl));

%% Step 2: Run MX-ONMTF (trace minimization, no ground truth)

results = MX_ONMTF(Al, k, kc, kpl, ...
    'runs', 20, ...
    'max_iter', 2000, ...
    'assign_perc', 0.8);

%% Step 3: Examine results

fprintf('\nCommunity assignments (first 10 nodes per layer):\n');
for l = 1:L
    fprintf('  Layer %d: %s\n', l, num2str(results.Clusters{1}((l-1)*n+1:(l-1)*n+10)'));
end

%% Step 4: Evaluate with Modularity Density (optional)

for l = 1:L
    [~,Il] = max(results.H_best{1}, [], 2);
    [md, md_norm] = ModDen(k(l), Al{l}, Il);
    fprintf('  Layer %d Modularity Density: %.4f (normalized: %.4f)\n', l, md, md_norm);
end

%% Helper function: planted partition model
function A = generate_planted_partition(labels, p_in, p_out)
    n = length(labels);
    A = zeros(n);
    for i = 1:n
        for j = i+1:n
            if labels(i) == labels(j)
                A(i,j) = rand() < p_in;
            else
                A(i,j) = rand() < p_out;
            end
        end
    end
    A = A + A';
end
