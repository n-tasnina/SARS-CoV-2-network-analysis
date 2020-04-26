import rpy2.robjects as ro
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

from rpy2.robjects.conversion import localconverter

clusterProfiler = importr('clusterProfiler')
annotationHub = importr('AnnotationHub')
orgdb = importr('org.Hs.eg.db')

#
# k= 200
# geneData = pd.read_csv('geneList',sep = '\t')
# geneList = geneData['prot'][0:]
# with localconverter(ro.default_converter + pandas2ri.converter):
#   geneList = ro.conversion.py2rpy(geneList)
#
# gene = geneData['prot'][0:k]
# with localconverter(ro.default_converter + pandas2ri.converter):
#   gene = ro.conversion.py2rpy(gene)
#
# print(type(gene))
# print(gene)
#
# ego = clusterProfiler.enrichGO(gene  = gene,
# 		        universe      =   geneList ,
#               	keyType       = 'UNIPROT',
#                 OrgDb         = orgdb,
#                 ont           = "BP",
#                 pAdjustMethod = "BH",
#                 pvalueCutoff  = 0.01,
#                 qvalueCutoff  = 0.05)

ro.r('''
        library('org.Hs.eg.db')
        k=200
        
        geneData = read.table('geneList',sep = '\t')
        #head(geneData)
        
        geneList <- geneData[,3]
        ##geneList
        
        
        
        names(geneList) <- as.character(geneData[,2])
        
        
        geneList <- sort(geneList, decreasing = TRUE)
        
        ##head(geneList)
        
        gene <- names(geneList)[1:k]
        
        
        ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                        keyType       = 'UNIPROT',
                        OrgDb         = org.Hs.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05)
        head(summary(ego))
        
        write.csv(ego, 'enrichGO.csv')
        ego <- setReadable(ego, OrgDb = org.Hs.eg.db)

        dotplot(ego, showCategory=30)
        
        
    ''')

#
# ego <- clusterProfiler.setReadable(ego, OrgDb = orgdb)
#
# clusterProfiler.dotplot(ego, showCategory=30)


