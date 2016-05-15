function [ Mean_MRR_result , MRR_result ] = MRR( Result , mmu_mgi_mp_wiped , mmu_mgi_mp , top_n )

test = mmu_mgi_mp - mmu_mgi_mp_wiped;

[m,~] = size(Result);
[~,index]=sort(Result,2,'descend');
MRR_result = zeros(m,1);

%iterate all genes
for i = 1:m    
    ranked = test(i,(index(i,:)));
%     index_2 = find(ranked == 2);
%     ranked(index_2)=[];
    
    if sum(rank(1:top_n)) ~=0
        tmp = find(ranked);
        MRR_result(i,1) = 1/tmp(1);
    end
    
end

Mean_MRR_result = sum(MRR_result)/m;


end

