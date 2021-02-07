library(frustratometeR)
OrderList <- c("1a0i.pdb")
Pdb_mut <- dir_frustration(PdbsDir = "/home/maria/Documentos/FrustraPocket/job.1a0i/",Chain = "A", OrderList = OrderList, Mode = "mutational", ResultsDir = "/home/maria/Documentos/FrustraPocket/job.1a0i/A/")