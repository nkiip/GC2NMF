path(path,'data')
path(path, '../evaluation')
%load('data670.mat', 'mmu_mp_mgi');%mouse
%mmu_mgi_mp = mmu_mp_mgi';%mouse
load('hsa_data3193.mat', 'hsa_hp_ncbi');%human
hsa_ncbi_hp = hsa_hp_ncbi';%human

test_set_percent = 0.2;
train_set_percent = [0.2,0.4,0.6,0.8,1];
max_iter = 10;
top_n_set = [200,400,600,800,1000];
auc_map_ndcg_top_n_set = zeros(length(train_set_percent),length(top_n_set),7);
tmp = zeros(1,length(top_n_set), 7);

for m = 1:length(train_set_percent)
    for t = 1:max_iter
        %* random choose test & validation & train data_set
        [train_vec, probe_vec, genes, mps, train_set] = data_trans(hsa_ncbi_hp, test_set_percent ,train_set_percent(m));

        %* run pmf;
        %[reslult, w1_P1, w1_M1] = pmf(train_vec, probe_vec, genes, mps);% for pmf

        %* run bpmf;
        [~, w1_P1, w1_M1] = pmf(train_vec, probe_vec, genes, mps);% 1-bpmf
        result = bayespmf(train_vec, probe_vec, w1_P1, w1_M1, genes, mps);% 2-bpmf

        %* evaluation
        for i = 1:length(top_n_set)
            mmu_mgi_mp_predict = result;
            mmu_mgi_mp_wiped = train_set;
            mmu_mgi_mp_test_set = hsa_ncbi_hp - train_set;
            [Mean_AUC_result , AUC_result] = AUC_new(mmu_mgi_mp_predict , mmu_mgi_mp_wiped , hsa_ncbi_hp , top_n_set(i));
            [map_value, map] = MAP(mmu_mgi_mp_predict, mmu_mgi_mp_test_set , top_n_set(i));
            [ndcg_value, ndcg] = NDCG(mmu_mgi_mp_predict, mmu_mgi_mp_test_set , top_n_set(i));
            [Mean_F1_result, f1 , Mean_Recall_result , Recall , Mean_Precision_result , Precision] = F_value(mmu_mgi_mp_predict , mmu_mgi_mp_wiped , hsa_ncbi_hp , top_n_set(i));
            [Mean_MRR_result, mrr] = MRR(mmu_mgi_mp_predict , mmu_mgi_mp_wiped , hsa_ncbi_hp , top_n_set(i));
            tmp(1,i,:) = [Mean_AUC_result , map_value, ndcg_value,Mean_F1_result,Mean_Recall_result,Mean_Precision_result,Mean_MRR_result];
        end
        
        auc_map_ndcg_top_n_set(m,:,:) = auc_map_ndcg_top_n_set(m,:,:) + tmp(1,:,:);
    end
end
auc_map_ndcg_top_n_set = auc_map_ndcg_top_n_set/max_iter;

directory='../4_result/GNMF/';
if(~exist(directory,'dir'))
 mkdir(directory);
end
%store the result
datetime=fix(clock);
t='';
for i=1:length(datetime)
    t=[t num2str(datetime(i))];
end
addresult = [directory 'result_bpmf_' t '.mat'];
save(addresult,'auc_map_ndcg_top_n_set','max_iter');
% save(result,'batch_folds','max_ites','epsilon','-append');
% save(result,'K','top_n_set','-append');
% save(result,'evaluate_result','obj_val_result','mean_evaluate_result','-append');
% save(result,'test_set_percent','train_set_percent','-append');
% save(result,'U_result','V_result','-append');