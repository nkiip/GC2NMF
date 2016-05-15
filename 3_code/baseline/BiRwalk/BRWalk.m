path(path, '../evaluation');
% load('data670');
% mmu_mgi_mp = mmu_mp_mgi'; 
load('hsa_data3193.mat')
hsa_ncbi_hp = hsa_hp_ncbi';

test_percent = 0.2;
hsa_ppi = full(hsa_ppi);
%%该函数用于把一个稀疏矩阵（sparse matrix）转换成一个全矩阵（full matrix）
ppi = diag(1/sum(hsa_ppi, 2)) * hsa_ppi;
pp = hsa_hp_hp * diag(1./sum(hsa_hp_hp, 1));
alpha = 1;
max_iter = 10;
top_n_set = [200,400,600,800,1000];
auc_map_ndcg_top_n_set = zeros(length(top_n_set),7);
tmp = zeros(length(top_n_set), 7);

for t = 1:max_iter
    %* random choose test & validation & train data_set
    test = Random_Choose_Test_Set(hsa_ncbi_hp, test_percent);
    %%设置测试集
    train = hsa_ncbi_hp - test;
    
    %* run BRwalk
    result = alpha * ppi * train  + (1-alpha) * train;
    
    %* evaluation
    for i = 1:length(top_n_set)
        mmu_mgi_mp_predict = result;
        mmu_mgi_mp_wiped = train;
        mmu_mgi_mp_test_set = test;
        [Mean_AUC_result , AUC_result] = AUC_new(mmu_mgi_mp_predict , mmu_mgi_mp_wiped , hsa_ncbi_hp , top_n_set(i));
        [map_value, map] = MAP(mmu_mgi_mp_predict, mmu_mgi_mp_test_set , top_n_set(i));
        [ndcg_value, ndcg] = NDCG(mmu_mgi_mp_predict, mmu_mgi_mp_test_set , top_n_set(i));
        [Mean_F1_result, f1 , Mean_Recall_result , Recall , Mean_Precision_result , Precision] = F_value(mmu_mgi_mp_predict , mmu_mgi_mp_wiped , hsa_ncbi_hp , top_n_set(i));
        [Mean_MRR_result, mrr] = MRR(mmu_mgi_mp_predict , mmu_mgi_mp_wiped , hsa_ncbi_hp , top_n_set(i));
        tmp(i,:) = [Mean_AUC_result , map_value, ndcg_value,Mean_F1_result,Mean_Recall_result,Mean_Precision_result,Mean_MRR_result];
    end
    auc_map_ndcg_top_n_set = auc_map_ndcg_top_n_set + tmp;
end
auc_map_ndcg_top_n_set = auc_map_ndcg_top_n_set/max_iter;
