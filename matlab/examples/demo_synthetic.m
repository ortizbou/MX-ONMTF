%% MX-ONMTF Demo: Synthetic Multiplex Network
% This script demonstrates how to run MX-ONMTF on synthetic multiplex
% networks with known ground truth community labels.
%
% Reference:
%   M. Ortiz-Bouza and S. Aviyente, "Community Detection in Multiplex
%   Networks Based on Orthogonal Nonnegative Matrix Tri-Factorization,"
%   IEEE Access, vol. 12, pp. 6423-6436, 2024.

%% Add paths
addpath('..');
addpath('../helpers');

%% Step 0: Generate a synthetic multiplex network
% In practice, replace this section with your own data loading.
% Here we create a simple 2-layer network with known community structure.

n = 60;          % number of nodes
L = 2;           % number of layers
kc_true = 2;     % true number of common communities

% Ground truth: 3 communities per layer, 2 common + 1 private each
GT1 = [ones(20,1); 2*ones(20,1); 3*ones(20,1)];
GT2 = [ones(20,1); 2*ones(20,1); 4*ones(20,1)];  % community 4 is private to layer 2

% Generate adjacency matrices from planted partition model
p_in = 0.4;   % intra-community edge probability
p_out = 0.05; % inter-community edge probability

A1 = generate_planted_partition(GT1, p_in, p_out);
A2 = generate_planted_partition(GT2, p_in, p_out);

% Package as cell arrays (single realization)
Al = {A1, A2};
GTl = {GT1, GT2};

%% Step 1: Find the number of communities (optional)
% If the number of communities is known, skip this step and set kc, k, kpl
% directly.

[kc, k, kpl] = findingk(Al);

fprintf('Detected communities:\n');
fprintf('  Common communities (kc): %d\n', kc);
fprintf('  Total per layer (k):     [%s]\n', num2str(k));
fprintf('  Private per layer (kpl): [%s]\n', num2str(kpl));

%% Step 2: Run MX-ONMTF with ground truth (NMI-based selection)

results = MX_ONMTF(Al, k, kc, kpl, ...
    'ground_truth', GTl, ...
    'realizations', 1, ...
    'runs', 20, ...
    'max_iter', 1000);

fprintf('\nResults:\n');
fprintf('  Best NMI: %.4f\n', results.averagedNMI);

%% Step 3: Examine community assignments

fprintf('\nCommunity assignments (first 10 nodes per layer):\n');
for l = 1:L
    fprintf('  Layer %d: %s\n', l, num2str(results.Clusters{1}((l-1)*n+1:(l-1)*n+10)'));
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
    A = A + A';  % symmetrize
end
