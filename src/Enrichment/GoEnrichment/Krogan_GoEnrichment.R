library('clusterProfiler')
library('AnnotationHub')
library('org.Hs.eg.db')


dirname = "/media/tassnina/Study/VT/Research_Group/SARSCOV2/SARS-CoV-2-network-analysis/fss_inputs/pos-neg/2020-03-sarscov2-human-ppi"
Krogan_protein_file="/media/tassnina/Study/VT/Research_Group/SARSCOV2/SARS-CoV-2-network-analysis/fss_inputs/pos-neg/2020-03-sarscov2-human-ppi/pos.txt"


Krogan_protein_set = read.table(Krogan_protein_file, sep = '\t')


Krogan_protein_names <- Krogan_protein_set[,1]


#####BIOLOGICAL PROCESS ENRICHMENT

ego_BP <- enrichGO(gene          = Krogan_protein_names,
	      	keyType       = 'UNIPROT',
	        OrgDb         = org.Hs.eg.db,
	        ont           = "BP",
	        pAdjustMethod = "BH",
	        pvalueCutoff  = 0.01,
	        qvalueCutoff  = 0.05)

output_file = paste(dirname, 'Krogan_enrichGO_BP.csv',sep='/')

write.table(ego_BP,output_file,sep='\t' )


##### CELLULAR COMPONENT ENRICHMENT
ego_CC <- enrichGO(gene          = Krogan_protein_names,
	      	keyType       = 'UNIPROT',
	        OrgDb         = org.Hs.eg.db,
	        ont           = "CC",
	        pAdjustMethod = "BH",
	        pvalueCutoff  = 0.01,
	        qvalueCutoff  = 0.05)

output_file = paste(dirname, 'Krogan_enrichGO_CC.csv',sep='/')

write.table(ego_CC,output_file, sep='\t')


##### MOLECULAR FUNCTION ENRICHMENT
ego_MF <- enrichGO(gene          = Krogan_protein_names,
	      	keyType       = 'UNIPROT',
	        OrgDb         = org.Hs.eg.db,
	        ont           = "MF",
	        pAdjustMethod = "BH",
	        pvalueCutoff  = 0.01,
	        qvalueCutoff  = 0.05)

output_file = paste(dirname, 'Krogan_enrichGO_MF.csv',sep='/')

write.table(ego_MF,output_file, sep='\t' )


