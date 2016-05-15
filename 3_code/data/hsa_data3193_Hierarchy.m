path(path,'../../2_useful_data')
clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-load-%%%%%%%%%%%%%%%%%%%%%%%%%%
%phenotype-gene association matrix
load('matrix.mat','hsa_hp_ncbi','hsa_ncbi_ncbi');
load('index.mat','hsa_hp_ncbi_row','hsa_hp_ncbi_col','hsa_ncbi_ncbi_index');
%PPI gene-gene relation matrix
%load('ppi.mat','ppi','mgi_id');
hsa_ppi_ncbi_id = hsa_ncbi_ncbi_index;

%Kegg gene pathways
load('pathway.mat','pathway_hsa_ncbi');
hsa_genes_pathways_origin = unique(pathway_hsa_ncbi);
hsa_pathway_ncbi_origin = pathway_hsa_ncbi;
%the relation between phenotype and phenotype
load('relation.mat','rank_hp_hp','hp_id');

%filter out genes, which both occur in gene-phenotype matrix and PPI
hsa_hp_hp = 1+log(rank_hp_hp);
hsa_hp_hp(hsa_hp_hp<0) = 0;

[C,ia,ib] = intersect(hsa_ppi_ncbi_id, hsa_hp_ncbi_col);
hsa_ncbi_ids = C;
hsa_ppi = hsa_ncbi_ncbi(ia,ia);
hsa_hp_ncbi_origin = hsa_hp_ncbi;
hsa_hp_ncbi = hsa_hp_ncbi(:,ib);

%filter out phenotypes, which both occur in gene-phenotype matrix and
%phenotype-phenotype matrix
[hsa_hp_ids,ia,ib] = intersect(hsa_hp_ncbi_row,hp_id);
hsa_hp_ncbi = hsa_hp_ncbi(ia,:);
hsa_hp_hp_origin = hsa_hp_hp;
hsa_hp_hp = hsa_hp_hp(ib,ib);
%filter the empty genes and phenotypes
hsa_hp_ncbi = hsa_hp_ncbi(sum(hsa_hp_ncbi,2)>0,sum(hsa_hp_ncbi,1)>0);
hsa_ncbi_ids = hsa_ncbi_ids(sum(hsa_hp_ncbi,1)>0);
hsa_hp_ids = hsa_hp_ids(sum(hsa_hp_ncbi,2)>0);
hsa_ppi = hsa_ppi(sum(hsa_hp_ncbi,1)>0,sum(hsa_hp_ncbi,1)>0);
hsa_hp_hp = hsa_hp_hp(sum(hsa_hp_ncbi,2)>0,sum(hsa_hp_ncbi,2)>0);
%using the new selected genes for each pathway
[rows,~] = size(hsa_pathway_ncbi_origin);
hsa_pathway_ncbi = zeros(rows,length(hsa_ncbi_ids));
for i = 1:rows
    [C2,ia2,~] = intersect(hsa_ncbi_ids, hsa_pathway_ncbi_origin(i,:));
    hsa_pathway_ncbi(i,ia2) = 1;
end
save('hsa_data3193_Hierarchy.mat','hsa_hp_ncbi');
save('hsa_data3193_Hierarchy.mat','hsa_ncbi_ids','hsa_hp_ids','-append');
save('hsa_data3193_Hierarchy.mat','hsa_pathway_ncbi','hsa_ppi','-append');
save('hsa_data3193_Hierarchy.mat','hsa_hp_hp','-append');



