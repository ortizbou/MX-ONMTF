%% Applying Symmetric Nonnegative Matrix 3-Factorization with orthogonal
%% restriction to each layer separately.
%%
%% When ground_truth is provided, uses NMI to select the best solution.
%% When omitted, uses Modularity Density.
function [U1_best, best_metric] = ONMTF(A, k, opts)

%  Input:
%         A:              adjacency matrix
%         k:              number of communities
%
%  Optional name-value arguments:
%         ground_truth:   ground truth labels. When provided, NMI is used
%                         for solution selection. When omitted, Modularity
%                         Density is used. Default: []
%         eta:            learning rate in (0,1], default 0.5
%         runs:           number of independent runs, default 50
%         max_iter:       maximum iterations per run, default 1000
%         conv_thresh:    convergence threshold, default 1e-3
%
%  Output:
%         U1_best:        best factor matrix across runs
%         best_metric:    best NMI or Modularity Density value
%
%  Author: Meiby Ortiz-Bouza
%  Address: Michigan State University, ECE
%  email: ortizbou@msu.edu

arguments
    A
    k
    opts.ground_truth       = []
    opts.eta          (1,1) double = 0.5
    opts.runs         (1,1) double = 50
    opts.max_iter     (1,1) double = 1000
    opts.conv_thresh  (1,1) double = 1e-3
end

eta = opts.eta;
use_nmi = ~isempty(opts.ground_truth);
[~,n] = size(A);

for j=1:opts.runs
    %% Initializing U1, C1
    U1 = rand([n,k]);
    C1 = diag(rand([k,1]));

    %% Iterative update
    for i=1:opts.max_iter
        U1new = U1.*((A*U1*C1')./(U1*U1'*A*U1*C1')).^eta;
        C1new = C1.*((U1'*A*U1)./(U1'*U1*C1*(U1'*U1))).^eta;
        if all(isnan(U1new),'all')==1
            U1new = rand([n,k]);
        end
        if all(isnan(C1new),'all')==1
            C1new = rand([k,k]);
        end
        if (all(norm(U1-U1new)<opts.conv_thresh,'all') && all(norm(C1-C1new)<opts.conv_thresh,'all'))
            break
        end
        U1 = U1new;
        C1 = C1new;
    end

    %% Finding Clusters and Computing metric
    [~,I1] = max(U1,[],2);

    if use_nmi
        %% NMI-based selection (synthetic)
        metric_j = getNMI(opts.ground_truth, I1);
        metric_all(j) = metric_j;
        max_metric = max(metric_all);

        if max_metric == metric_j
            U1_best = U1;
        end

        if abs(max_metric - 1) < 1e-3
            break
        end
    else
        %% Modularity Density selection (real)
        metric_j = ModDen(k, A, I1);
        metric_all(j) = metric_j;
        max_metric = max(metric_all);

        if max_metric == metric_j
            U1_best = U1;
        end
    end
end

best_metric = max(metric_all);

end
