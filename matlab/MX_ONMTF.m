%%%%% Multiplex community detection with common and private communities
%%%%% L layers
function results = MX_ONMTF(Alr, k, kc, kpl, opts)
%  Multiplex Orthogonal Nonnegative Matrix Tri-Factorization.
%
%  When ground_truth is provided, uses NMI to select the best solution
%  (synthetic mode). When omitted, uses trace minimization (real mode).
%
%  Input:
%         Alr:  cell array of adjacency matrices per layer (single realization),
%               or cell array of cell arrays (multiple realizations)
%         k:    L vector with the total number of communities per layer
%         kc:   number of common communities
%         kpl:  L vector with the number of private communities per layer
%
%  Optional name-value arguments:
%         ground_truth:   cell array of ground truth labels per layer (single realization),
%                         or cell array of cell arrays (multiple realizations).
%                         When provided, NMI is used for solution selection.
%                         When omitted, trace minimization is used. Default: {}
%         eta:            learning rate in (0,1], default 0.5
%         runs:           number of independent runs per realization, default 20
%         realizations:   number of network realizations, default 1
%         max_iter:       maximum iterations per run, default 1000
%         conv_thresh:    convergence threshold, default 1e-3
%         assign_mode:    community assignment mode ('flex' or 'same'), default 'same'
%         assign_perc:    threshold for common community assignment, default 0.8
%
%  Output: struct with fields:
%         H_best:           best common community matrix per realization
%         Hl_best:          best private community matrices per realization
%         Clusters:         community labels
%         NMI_realization:  (NMI mode only) best NMI per realization
%         averagedNMI:      (NMI mode only) mean NMI across realizations
%         stdNMI:           (NMI mode only) std NMI across realizations
%
%  Author: Meiby Ortiz-Bouza
%  Address: Michigan State University, ECE
%  email: ortizbou@msu.edu
%
%  Reference:
%  [1] "Community detection in multiplex networks based on orthogonal
%       nonnegative matrix tri-factorization"
%       Authors: Meiby Ortiz-Bouza and Selin Aviyente

arguments
    Alr
    k
    kc
    kpl
    opts.ground_truth       cell   = {}
    opts.eta          (1,1) double = 0.5
    opts.runs         (1,1) double = 20
    opts.realizations (1,1) double = 1
    opts.max_iter     (1,1) double = 1000
    opts.conv_thresh  (1,1) double = 1e-3
    opts.assign_mode  string       = "same"
    opts.assign_perc  (1,1) double = 0.8
end

eta = opts.eta;
use_nmi = ~isempty(opts.ground_truth);

for r=1:opts.realizations
    %% Extract layers for this realization
    if opts.realizations > 1
        Al = Alr{r};
    else
        Al = Alr;
    end
    L = size(Al,2);
    n = size(Al{1},2);

    if use_nmi
        if opts.realizations > 1
            GTl = opts.ground_truth{r};
        else
            GTl = opts.ground_truth;
        end
        GTSup = vertcat(GTl{:});
    end

    %% Running the code multiple times
    for j=1:opts.runs
        %% Initializing Hl, H, Sl, Gl
        for l=1:L
            Hl{l}=rand([n,kpl(l)]);
            Sl{l}=diag(rand(kc,1));
            Gl{l}=diag(rand(kpl(l),1));
        end
        H=rand([n,kc]);

        %% Iterative update
        for i=1:opts.max_iter
            for l=1:L
                Hlnew{l}=Hl{l}.*((Al{l}*Hl{l}*Gl{l} + Hl{l}*Hl{l}'*H*Sl{l}*H'*Hl{l}*Gl{l}')./(H*Sl{l}*H'*Hl{l}*Gl{l} + Hl{l}*Hl{l}'*Al{l}*Hl{l}*Gl{l})).^eta;
                if all(isnan(Hlnew{l}),'all')==1
                    Hlnew{l}=rand([n,kpl(l)]);
                end
            end

            numH=0;
            denH=0;
            for l=1:L
                numH=numH + Al{l}*H*Sl{l} + H*H'*Hl{l}*Gl{l}'*Hl{l}'*H*Sl{l};
                denH=denH + Hl{l}*Gl{l}'*Hl{l}'*H*Sl{l} +  H*H'*Al{l}*H*Sl{l};
            end
            Hnew=H.*(numH./denH).^eta;
            if all(isnan(Hnew),'all')==1
                Hnew=rand([n,kc]);
            end

            for l=1:L
                Slnew{l}=Sl{l}.*((H'*Al{l}*H)./(H'*H*Sl{l}*(H'*H) + H'*Hl{l}*Gl{l}*Hl{l}'*H)).^eta;
                if all(isnan(Slnew{l}),'all')==1
                    Slnew{l}=diag(rand(kc,1));
                end

                Glnew{l}=Gl{l}.*((Hl{l}'*Al{l}*Hl{l})./(Hl{l}'*Hl{l}*Gl{l}*(Hl{l}'*Hl{l}) + Hl{l}'*H*Sl{l}*H'*Hl{l})).^eta;
                if all(isnan(Glnew{l}),'all')==1
                    Glnew{l}=diag(rand([kpl(l),1]));
                end
            end

            Hres=cellfun(@minus,Hl,Hlnew,'Un',0);
            Sres=cellfun(@minus,Sl,Slnew,'Un',0);
            if (all(norm(Hres{l})<opts.conv_thresh,'all') && ...
                 all(norm(Sres{l})<opts.conv_thresh,'all') ...
            && all(norm(H-Hnew)<opts.conv_thresh,'all'))
                break
            end

            Hl=Hlnew;
            H=Hnew;
            Sl=Slnew;
            Gl=Glnew;
        end

        %% Finding Clusters
        [Il,ClustersSupra] = assigncomm(Al,H,Hl,kc,k,kpl,L,opts.assign_mode,opts.assign_perc);

        %% Solution selection
        if use_nmi
            %% NMI-based selection (synthetic)
            for l=1:L
                NMIlj{l}=getNMI(GTl{l},Il{l});
            end
            NMIsupj=getNMI(GTSup,ClustersSupra);
            NMI_sup(j)=NMIsupj;
            maxNMI=max(NMI_sup);

            if maxNMI==NMIsupj
                results.H_best{r}=H;
                results.Hl_best{r}=Hl;
                results.Clusters=ClustersSupra;
            end

            if norm(maxNMI-1)<1e-3
                break
            end
        else
            %% Trace minimization (real)
            for l=1:L
                tr(l)=norm(Al{l}-(H*Sl{l}*H'+Hl{l}*Gl{l}*Hl{l}'),'fro')+trace(Hl{l}'*Hl{l}-eye(kpl(l)));
            end
            trall(j)=sum(tr)+trace(H'*H-eye(kc));
            mintrace=min(trall);

            if mintrace==trall(j)
                results.H_best{r}=H;
                results.Hl_best{r}=Hl;
                results.Clusters{r}=ClustersSupra;
            end
        end
    end

    if use_nmi
        [results.NMI_realization(r),~]= max(NMI_sup);
        clear NMI_sup
    end
end

if use_nmi
    results.averagedNMI = mean(results.NMI_realization);
    results.stdNMI = std(results.NMI_realization);
end

end
