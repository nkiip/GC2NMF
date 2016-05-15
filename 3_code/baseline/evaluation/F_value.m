function [ Mean_F1_result , F1_result , Mean_Recall_result,Recall_result, Mean_Precision_result , Precision_result] = F_value( Result , mmu_mgi_mp_wiped , mmu_mgi_mp , top_n )

test = mmu_mgi_mp - mmu_mgi_mp_wiped;

[m,~] = size(Result);
[~,index]=sort(Result,2,'descend');
F1_result = zeros(m,1);
Recall_result = zeros(m,1);
Precision_result = zeros(m,1);

for i = 1:m    
    ranked = test(i,(index(i,:)));
%   index_2 = find(ranked == 2);
%   ranked(index_2)=[];
    Recall_result(i,1) = sum(ranked(1:top_n))/sum(test(i,:));
    Precision_result(i,1) = sum(ranked(1:top_n))/top_n;
    if Recall_result(i,1)+Precision_result(i,1)~=0
    F1_result(i,1) = Recall_result(i,1)*Precision_result(i,1)/(Recall_result(i,1)+Precision_result(i,1));
    end
end

Mean_F1_result = sum(F1_result)/m;
Mean_Recall_result = sum(Recall_result)/m;
Mean_Precision_result = sum(Precision_result)/m;

end

