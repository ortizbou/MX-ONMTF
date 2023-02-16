%%%%% Multiplex with common communities
%%%%% L layers
function [H_bestNMI,Hl_bestNMI]=MX_ONMTF_real(Alr,k,kc,kpl)
% Alr is cell array with R realizations of the networks. Each cell contains another L cells array with the adjacency matrices of each L layer
% k is an L vector with the total number of communities per layer
% kc is the number of common communities
% kpl is an L vector with the number of private communities per layer
% communities membership matrices

eta=1; %% learning rate between (0,1]
realizations=1;
% n=256;

for r=1:realizations
Al=Alr;
L=size(Al,2);
n=size(Al{1},2); 

%% Running the code multiple times and finding NMI
runs=20;

for j=1:runs
    %% Initializing H1,H2,H,S1,S2,G1,G2
    for l=1:L
    Hl{l}=rand([n,kpl(l)]);
    Sl{l}=diag(rand(kc,1));
    Gl{l}=diag(rand(kpl(l),1));
    end
    H=rand([n,kc]);

    %% Iterative update
    for i=1:1000
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
        if (all(norm(Hres{l})<1e-4,'all') && ...
             all(norm(Sres{l})<1e-4,'all') ...
        && all(norm(H-Hnew)<1e-4,'all'))
            break
        end    
    
        Hl=Hlnew;      
        H=Hnew;
        Sl=Slnew;
        Gl=Glnew;
        it=i;
    end

%% Finding Clusters
% 
% for l=1:L
%     Hl{l}=[H,Hl{l}];
%     for m=1:(kc-(k(l)-kpl(l)))
%     pos=patchmult(Al{l},Hl{l},kc);
%     Hl{l}(:,pos(m))=0;
%     end
%     [~,Il{l}] = max(Hl{l},[],2);
% end
% 
% for l=1:L
%     for m=1:kpl(l)
%     a=find(Il{l}==(kc+m));
%     Il{l}(a)=Il{l}(a)+sum(k(1:(l-1)));
%     end
% end
% 
% ClustersSupra=[vertcat(Il{:})]; 
for l=1:20
    Il{l}=zeros(n,1);
end
Hall=horzcat(Hl{:});
for g=1:kc
    for nodes=1:n
    if all(H(nodes,g)>Hall(nodes,:))
    for l=1:L    
        Il{l}(nodes)=g;
    end
    end
    end
end
nodes=find(Il{1}==0);
    for l=1:L
        il=zeros(n,1);
        [~,il(nodes)] = max(Hl{l}(nodes,:),[],2);
        Il{l}(nodes)= il(nodes)+kc+sum(k(1:(l-1)));
    end  
ClustersSupra=[vertcat(Il{:})]'; 

%% Modularity
for l=1:L
[ModDenj{l},ModDenNormj{l}]=ModDen(k(l),Al{l},Il{l});
end
Asup=zeros(L*n);
for l=1:L
    Asupra(1+n*(l-1):n*l,1+n*(l-1):n*l)=Al{l};
end
[ModDenSupj,ModDenNormSupj]=ModDen(sum(kpl)+kc,Asupra,ClustersSupra);
Mod_DenSup(j)=ModDenSupj;
maxModDenSup=max(Mod_DenSup);

if maxModDenSup==ModDenSupj
   H_bestNMI{r}=H; 
   Hl_bestNMI{r}=Hl;   
%    Sl_bestNMI{r}=Sl; 
%    Gl_bestNMI{r}=Gl; 
end   

end

[ModDen_realization(r),~]= max(Mod_DenSup);

end


end
