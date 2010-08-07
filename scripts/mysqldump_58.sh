mysqldump -uensadmin -pensembl -hens-livemirror -P3306 --opt --routines ensembl_compara_58 meta ncbi_taxa_node ncbi_taxa_name genome_db sequence sequence_cds member subset subset_member protein_tree_node protein_tree_member protein_tree_tag protein_tree_stable_id species_set species_set_tag method_link method_link_species_set > ensembl_compara_58.sql

mysqldump -uensadmin -pensembl -hens-livemirror -P3306 --opt --routines ensembl_compara_58 conservation_score constrained_element dnafrag > ensembl_compara_58.sql
