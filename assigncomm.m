function [Il,ClustersSupra]=assigncomm(Al,H,Hl,kc,k,kpl,L,mode)

[n,~]=size(H);

for l=1:L
    Hc{l}=H;
    m=(kc-(k(l)-kpl(l)));
    pos=patchmult(Al{l},H,kc);
    Hc{l}(:,pos(1:m))=0;
    Il{l}=zeros(n,1);
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

    for nodes=1:n
        for g=1:kc
        for l=1:L
            H0l(:,l)=Hc{l}(:,g);
        end
        layers=find(sum(H0l)~=0);
        Hall=horzcat(Hl{layers});
            if (all(H(nodes,g)>Hall(nodes,:)) && all(H(nodes,g)>=H(nodes,:)))
                for l=layers 
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