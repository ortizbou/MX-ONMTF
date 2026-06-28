function N = rnorm(M)
%% normalize the rows of matrix M
N = M;

for i = 1:size(M,1)
    if norm(M(i,:)) ~= 0
        N(i,:) = M(i,:)./norm(M(i,:));
    end
end