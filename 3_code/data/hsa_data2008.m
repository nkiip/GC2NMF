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
load('relation.mat','hsa_hp_hp','hp_id');

[C,~,ib] = intersect(hsa_pathway_ncbi_origin, hsa_hp_ncbi_col);
hsa_ncbi_ids = C;
hsa_hp_ncbi_origin = hsa_hp_ncbi;
hsa_hp_ncbi = hsa_hp_ncbi(:,ib);

%filter phenotype
[hsa_hp_ids,ia,ib] = intersect(hsa_hp_ncbi_row,hp_id);
hsa_hp_ncbi = hsa_hp_ncbi(ia,:);
hsa_hp_hp_origin = hsa_hp_hp;
hsa_hp_hp = hsa_hp_hp_origin(ib,ib);
%using the new selected genes for each pathway
[rows,~] = size(hsa_pathway_ncbi_origin);
hsa_pathway_ncbi = zeros(rows,length(hsa_ncbi_ids));
for i = 1:rows
    [~,ia,~] = intersect(hsa_ncbi_ids, hsa_pathway_ncbi_origin(i,:));
    hsa_pathway_ncbi(i,ia) = 1;
end
%using the new selected genes for gene PPI
hsa_ppi = zeros(length(hsa_ncbi_ids));
[~,ia,ib] = intersect(hsa_ncbi_ids, hsa_ppi_ncbi_id);
hsa_ppi(ia,ia) = hsa_ncbi_ncbi(ib,ib);


save('hsa_data2008.mat','hsa_hp_ncbi');
save('hsa_data2008.mat','hsa_ncbi_ids','hsa_hp_ids','-append');
save('hsa_data2008.mat','hsa_pathway_ncbi','hsa_ppi','-append');
save('hsa_data2008.mat','hsa_hp_hp','-append');