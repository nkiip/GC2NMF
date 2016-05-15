path(path, '../../evaluation');
%load('data670.mat', 'mmu_mp_mgi');%mouse
%mmu_mgi_mp = mmu_mp_mgi';%mouse
load('hsa_data3193', 'hsa_hp_ncbi');
hsa_ncbi_hp = hsa_hp_ncbi';

test_set_percent = 0.2;
train_set_percent=[0.2 0.4 0.6 0.8 1];
validaton_set_percent = [0.8 0.6 0.4 0.2 0];
max_iter = 10;
inittype = 2;
k = 5;
top_n_set = [200,400,600,800,1000];
auc_map_ndcg_top_n_set = zeros(length(train_set_percent),length(top_n_set),7);
tmp = zeros(1,length(top_n_set), 7);


for m=1:length(train_set_percent)
    for t = 1:max_iter
        %* random choose test & validation & train data_set
        R_test = Random_Choose_Test_Set(hsa_ncbi_hp, test_set_percent);
        R_v =  Random_Choose_Test_Set(hsa_ncbi_hp - R_test, validaton_set_percent(m));
        R_train = hsa_ncbi_hp-R_test-R_v;



        %* run ppmf;
        result = ppmf(R_train, R_v, R_test, R_train, R_v, R_test, inittype, k);

        %* evaluation
        for i = 1:length(top_n_set)
            mmu_mgi_mp_predict = result;
            mmu_mgi_mp_wiped = hsa_ncbi_hp - R_test;
            mmu_mgi_mp_test_set = R_test;
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

directory='../5_result/GNMF/second';
if(~exist(directory,'dir'))
 mkdir(directory);
end
%store the result
datetime=fix(clock);
t='';
for i=1:length(datetime)
    t=[t num2str(datetime(i))];
end
addresult = [directory 'result_ppmf_' t '.mat'];
save(addresult,'auc_map_ndcg_top_n_set','max_iter');