function [Il,ClustersSupra]=assigncomm(Alr,H,Hl,kc,k,kpl,L,mode)

% [n,~]=size(H);
% H=H/norm(H);
% for l=1:L
% Hl{l}=H{l}/norm(Hl{l});
% end
Al=Alr;
n=size(Al{1},2); 

for l=1:L
    Hc{l}=H;
for m=1:(kc-(k(l)-kpl(l)))
    pos=patchmult(Al{l},H,kc);
    Hc{l}(:,pos(m))=0;
end
end


if strcmp(mode,'flex')
    for l=1:L
        Hlnew{l}=[Hc{l},Hl{l}];
        [~,Il{l}] = max(Hlnew{l},[],2);
    end    
    for l=1:L
        for m=1:kpl(l)
        a=find(Il{l}==(kc+m));
        Il{l}(a)=Il{l}(a)+sum(k(1:(l-1)));
        end
    end 
    ClustersSupra=[vertcat(Il{:})];

elseif strcmp(mode,'same')
    for l=1:L
        Il{l}=zeros(n,1);
    end
    Hall=horzcat(Hl{:});
    for g=1:kc
        for nodes=1:n
            for l=1:L 
%                 if all(Hc{l}(nodes,g)>Hall(nodes,:))
                if nnz(Hc{l}(nodes,g)>Hall(nodes,:))>0.8*size(Hall,2)
                Il{l}(nodes)=g;
                end
            end
        end
    end
    for l=1:L
    nodes=find(Il{l}==0);
         il=zeros(n,1);  
            if ~isempty(Hl{l})                       
                [~,il(nodes)] = max(Hl{l}(nodes,:),[],2);
                Il{l}(nodes)= il(nodes)+kc+sum(k(1:(l-1)));               
            end
    end
    ClustersSupra=[vertcat(Il{:})];   
else 
    error('error in labelcomm.m : mode must be ''flex'' or ''same''');
end


end