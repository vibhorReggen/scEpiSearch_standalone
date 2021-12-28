library(UniPath)

args = commandArgs(trailingOnly=TRUE)
print(args[1])
if(length(args) ==0) {
 stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
peakfile = args[1]
count_matrix = args[2]
count_mat = read.csv(count_matrix, header = FALSE, sep = ",")
write.table(count_mat,file="./mat_unipath.txt",sep=",", col.names = F, row.names = T, quote = F)
type = args[3]

if (type==1){
  globalaccess = global_access(peakfile,"./storage/scepisearch/HumanFiles/reference_peaks_human.txt","./storage/scepisearch/HumanFiles/global_accessibility_score_human.csv")
  #imputed_count = drimpute("count_matrix.csv")
  cmd = nearest_gene("./storage/scepisearch/HumanFiles/nearestGenes.pl",peakfile,"./storage/scepisearch/HumanFiles/refseq-hg19.txt","foreground_output")
  ##system call
  system(cmd)
  data("c2.cp.v6.1.symbols")
  scores = runGO(c2.cp.v6.1.symbols,"./storage/scepisearch/HumanFiles/Background_global_human","./mat_unipath.txt",method=1,"globalaccess","foreground_output",promoters = FALSE,dist=1000000,threshold=1.25)
  adjpva = -log2(scores$hypergeometeric+.001)
  
  }else {
    globalaccess = global_access(peakfile,"./storage/scepisearch/MouseFiles/reference_peak_mouse","./storage/scepisearch/MouseFiles/global_accessibility_score_mouse.csv")
    #imputed_count = drimpute("count_matrix.csv")
    cmd = nearest_gene("./storage/scepisearch/HumanFiles/nearestGenes.pl",peakfile,"./storage/scepisearch/MouseFiles/mm9_refGene.txt","foreground_output")
    ##system call
    system(cmd)
    data("c2.cp.v6.1.symbols")
    scores = runGO(c2.cp.v6.1.symbols,"./storage/scepisearch/MouseFiles/background_mouse","./mat_unipath.txt",method=1,"globalaccess","foreground_output",promoters = FALSE,dist=1000000,threshold=1.25)
    adjpva = -log2(scores$hypergeometeric+.001)
  }

  write.table(adjpva,file="./pathway_scores.txt",sep=" ", col.names = T, row.names = T, quote = F)
 

