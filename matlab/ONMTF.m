%% Applying Symmetric Nonnegative Matrix 3-Factorization with orthogonal restriction to each layer separately
%%%% 50 runs are performed and the U1 with the best NMI is chosen. This is
%%%% repeated over 100 realizations of the same network and the averaged
%%%% NMI is reported. Try different values of k
function [U1_bestNMI,NMI]=ONMTF(A,GT,k,opts)

%  Optional name-value arguments:
%         eta:            learning rate in (0,1], default 0.5
%         runs:           number of independent runs, default 50
%         max_iter:       maximum iterations per run, default 1000
%         conv_thresh:    convergence threshold, default 1e-3
%
%  Author: Meiby Ortiz-Bouza
%  Address: Michigan State University, ECE
%  email: ortizbou@msu.edu

arguments
    A
    GT
    k
    opts.eta          (1,1) double = 0.5
    opts.runs         (1,1) double = 50
    opts.max_iter     (1,1) double = 1000
    opts.conv_thresh  (1,1) double = 1e-3
end

eta = opts.eta;
[~,n]=size(A);

%% Running the code multiple times and finding NMI

for j=1:opts.runs
    %% Initializing U1,C1
U1=rand([n,k]);
C1=diag(rand([k,1]));

%% Iterative update
for i=1:opts.max_iter
    U1new=U1.*((A*U1*C1')./(U1*U1'*A*U1*C1')).^eta;
    C1new=C1.*((U1'*A*U1)./(U1'*U1*C1*(U1'*U1))).^eta;
    if all(isnan(U1new),'all')==1
        U1new=rand([n,k]);
    end
    if all(isnan(C1new),'all')==1
        C1new=rand([k,k]);
    end
    if (all(norm(U1-U1new)<opts.conv_thresh,'all') && all(norm(C1-C1new)<opts.conv_thresh,'all'))
        break
    end
    U1=U1new;
    C1=C1new;
end

%% Finding Clusters and Computing NMI
[~,I1] = max(U1,[],2);
NMI1j=getNMI(GT,I1);
NMI(j)=NMI1j;
maxNMI=max(NMI);


if maxNMI==NMI1j
   U1_bestNMI=U1;
end


if abs(maxNMI-1)<1e-3
    break
end


end

[NMI,~]= max(maxNMI);


end
