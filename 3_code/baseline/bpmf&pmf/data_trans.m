function [train_vec, probe_vec, genes, mps, train_set] = data_trans(matrix, test_set_percent ,train_set_percent) 
test_set = Random_Choose_Test_Set(matrix , test_set_percent);%random choose test_test;
train_set1 = matrix - test_set;%get train set;
train_set = Random_Choose_Test_Set(train_set1 , train_set_percent);
[genes, mps] = size(matrix);%gete genes and phenos indices involed in process;
train_vec = zeros(nnz(train_set), 3);% define a triple vector store train_set data;
probe_vec = zeros(nnz(test_set), 3);%define a triple vector store test_set data;
[train_vec(:, 1), train_vec(:, 2), train_vec(:, 3)] = find(train_set > 0);%data trains process;
[probe_vec(:, 1), probe_vec(:, 2), probe_vec(:, 3)] = find(test_set > 0);%data trains process;
%save('mmudata.mat', 'train_vec', 'probe_vec', 'genes', 'mps', 'test_set');%save;
