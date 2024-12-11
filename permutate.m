function K = permutate(C)

c = numel(C);
c_aug = C + 1;

K = zeros(prod(c_aug),numel(c_aug));

for i = 1:c
    
    G = [];
    
    if i~=c
        L1 = prod(c_aug(i+1:end));
    else
        L1 = 1;
    end
    
    if i~=1
        L2 = prod(c_aug(1:i-1));
    else
        L2 = 1;
    end
    
    for j = 0:C(i)
        G = [G;j*ones(L1,1)];
    end
    
    K(:,i) = repmat(G,L2,1);
end

% K(sum(K,2)>Res,:)=[];