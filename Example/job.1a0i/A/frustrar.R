library(frustratometeR)
OrderList <- c("1a0i.pdb")
Pdb_sr <- dir_frustration(PdbsDir = "/home/camila/Documentos/2Frustra/job.1a0i/",Chain = "A", OrderList = OrderList, Mode = "singleresidue", ResultsDir = "/home/camila/Documentos/2Frustra/job.1a0i/A/")