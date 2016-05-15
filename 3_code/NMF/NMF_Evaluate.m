function [auc_map_ndcg_top_n_set] = NMF_Evaluate(mmu_mgi_mp_test_set, mmu_mgi_mp, U, V, top_n_set)
%GROUP_NMF_EVALUATE Summary of this function goes here
%   Detailed explanation goes here
    auc_map_ndcg_top_n_set = zeros(length(top_n_set),7);
    for i = 1:length(top_n_set)
        mmu_mgi_mp_predict = U*V;
        mmu_mgi_mp_wiped = mmu_mgi_mp - mmu_mgi_mp_test_set;
        [Mean_AUC_result , AUC_result] = AUC_new(mmu_mgi_mp_predict , mmu_mgi_mp_wiped , mmu_mgi_mp , top_n_set(i));
        [map_value, map] = MAP(mmu_mgi_mp_predict, mmu_mgi_mp_test_set , top_n_set(i));
        [ndcg_value, ndcg] = NDCG(mmu_mgi_mp_predict, mmu_mgi_mp_test_set , top_n_set(i));
        [Mean_F1_result, f1 , Mean_Recall_result , Recall , Mean_Precision_result , Precision] = F_value(mmu_mgi_mp_predict , mmu_mgi_mp_wiped , mmu_mgi_mp , top_n_set(i));
        [Mean_MRR_result, mrr] = MRR(mmu_mgi_mp_predict , mmu_mgi_mp_wiped , mmu_mgi_mp , top_n_set(i));
        auc_map_ndcg_top_n_set(i,:) = [Mean_AUC_result , map_value, ndcg_value,Mean_F1_result,Mean_Recall_result,Mean_Precision_result,Mean_MRR_result];
    end

end

