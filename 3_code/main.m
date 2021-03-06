path(path,'evaluation')
path(path,'data')
clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-load-%%%%%%%%%%%%%%%%%%%%%%%%%%
load('data670.mat');
% load('data670.mat','mmu_mp_mgi');
% load('data670.mat','mmu_mgi_ids','mmu_mp_ids');
% load('data670.mat','mmu_pathway_mgi','mmu_ppi');
% load('data670.mat','mmu_mp_mp');

mmu_mgi_mp = mmu_mp_mgi';
experiment_times = 10;
batch_folds = 1;
max_ites = 10*batch_folds; 
K = 10;
epsilon = 0.1;
 %we use variable 'test_set_percent' to set the test set percentage
test_set_percent = 0.2;
%'train_set_percent' doesn't equal to '1-test_set_percent', it
%indicates the percentage of origianl training set to train
train_set_percent = [0.2,0.4,0.6,0.8,1];

lambda0 =[0,0.001,0.01,0.1,1,10];
lambda1 = [0,0.001,0.01,0.1,1,10];
lambda2 = [1];

top_n_set = [200,400,600,800,1000];
evaluate_result = cell(length(train_set_percent),experiment_times, length(lambda0), length(lambda1), length(lambda2));
obj_val_result = cell(length(train_set_percent),experiment_times, length(lambda0), length(lambda1), length(lambda2));
U_result =  cell(length(train_set_percent),experiment_times, length(lambda0), length(lambda1), length(lambda2));
V_result =  cell(length(train_set_percent),experiment_times, length(lambda0), length(lambda1), length(lambda2));
for m = 1:length(train_set_percent)
    for i = 1:length(lambda0)
        for j = 1:length(lambda1)
            for k = 1:length(lambda2)
                %for a set of parameter combination, we repeat the exp
                tmp = cell(experiment_times,1);
                tmp2 = cell(experiment_times,1);
                tmp3 = cell(experiment_times,1);
                tmp4 = cell(experiment_times,1);
                for ite = 1:experiment_times

                    [mmu_mgi_mp_test_set] = Random_Choose_Test_Set(mmu_mgi_mp , test_set_percent);
                    [mmu_mgi_mp_train_set] = Train_set(mmu_mgi_mp, mmu_mgi_mp_test_set, train_set_percent(m));

                    %train the model
                    [U, V, tmp2{ite,1}] = Group_NMF_Train(mmu_mgi_mp_train_set, mmu_pathway_mgi, mmu_ppi, mmu_mp_mp, ...,
                        lambda0(i), lambda1(j), lambda2(k), K, max_ites, epsilon, batch_folds);
                    %evaluate the model
                    [tmp{ite,1}] = Group_NMF_Evaluate(mmu_mgi_mp_test_set, mmu_mgi_mp, U, V, top_n_set);
                    tmp3{ite,1} = U;
                    tmp4{ite,1} = V;
                end 
                evaluate_result(m,:,i,j,k) = tmp;
                obj_val_result(m,:,i,j,k) = tmp2;
                U_result(m,:,i,j,k) = tmp3;
                V_result(m,:,i,j,k) = tmp4;

                m
                i
                j
                k
                datestr(now) 
            end
        end
    end
end
mean_evaluate_result = cell(length(train_set_percent),length(lambda0), length(lambda1), length(lambda2));
for m = 1:length(train_set_percent)
    for i = 1:length(lambda0)
        for j = 1:length(lambda1)
            for k = 1:length(lambda2)
                tmp = zeros(size(evaluate_result{1,1,1,1,1}));
                for ite = 1:experiment_times
                    tmp = tmp + evaluate_result{m,ite,i,j,k};
                end
                mean_evaluate_result{m,i,j,k} = tmp/experiment_times;
            end
        end
    end
end
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
result = [directory 'result_' t '.mat'];
save(result,'lambda0','lambda1','lambda2','experiment_times');
save(result,'batch_folds','max_ites','epsilon','-append');
save(result,'K','top_n_set','-append');
save(result,'evaluate_result','obj_val_result','mean_evaluate_result','-append');
save(result,'test_set_percent','train_set_percent','-append');
save(result,'U_result','V_result','-append');