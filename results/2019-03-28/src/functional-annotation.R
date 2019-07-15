rm(list = ls())

## FUNCTIONS ==================================================================
annotate_culex <- function(df){
  annotation <- 
    t(sapply(df$atrribute, function(x){
      # Extract basic information from df
      atribute <- unlist(strsplit(as.character(x), ";"))
      xref <- atribute[startsWith(atribute, "Dbxref")]
      ncbi_id <- unlist(strsplit(xref, ","))[2]
      ncbi_id <- unlist(strsplit(ncbi_id, ":"))[2]
      vectorbase_id <- unlist(strsplit(xref, ","))[1]
      vectorbase_id <- unlist(strsplit(vectorbase_id, ":"))[2]
      
      # Donwload the gene description from the NCBI
      description <- rentrez::entrez_fetch(db="gene", id = ncbi_id,
                                           rettype = "text")
      description <- unlist(lapply(strsplit(description, "\n"), 
                                   function(x) x[3]))
      
      # Retrieve the ID from OrthoDB
      orthodb_info <- 
        try(read.table(paste0("https://www.orthodb.org/tab?query=", 
                              ncbi_id, "&ncbi=1&species=7157&level=7157"), 
                       sep = "\t", header = T))
      if(any(class(orthodb_info) == "try-error", nrow(orthodb_info) == 0)){
        orthodb_id <- ""
      } else {
        orthodb_id <- unique(as.character(orthodb_info$pub_og_id))
      }
      
      # Retrieve information from UniprotKB
      uniprot_info <- 
        try(read.table(paste0("https://www.uniprot.org/uniprot/?query=", 
                              ncbi_id, "&format=tab"), 
                       sep = "\t", header = T), silent = T)
      if(any(class(uniprot_info) == "try-error", nrow(uniprot_info) == 0)){
        uniprot_id <- ""
        gos_id <- ""
        gos_description <- ""
        interpro_id <- ""
        interpro_description <- ""
      } else {
        ## Uniprot ID
        uniprot_id <- as.character(uniprot_info[1, 1])
        uniprot_info <- read.table(paste0("https://www.uniprot.org/uniprot/", 
                                          uniprot_id, ".txt"),
                                   header = F, sep = "\t", 
                                   stringsAsFactors = F, fill = T)
        ## GO terms
        gos <- uniprot_info[startsWith(uniprot_info[, 1], "DR   GO"),]
        gos_id <- unlist(lapply(strsplit(gos, "; "), function(x) x[2]))
        gos_description <- unlist(lapply(strsplit(gos, "; "), function(x) x[3]))
        if(length(gos) > 1){
          gos_id <- paste(gos_id, collapse = ", ")
          gos_description <- paste(gos_description, collapse = "; ")
        } else if(length(gos) == 0){
          gos_id <- ""
          gos_description <- ""
        }
        ## InterPro annotation
        interpro <- uniprot_info[startsWith(uniprot_info[, 1], "DR   InterPro"),]
        interpro_id <- unlist(lapply(strsplit(interpro, "; "), function(x) x[2]))
        interpro_description <- unlist(lapply(strsplit(interpro, "; "), function(x) x[3]))
        if(length(interpro) > 1){
          interpro_id <- paste(interpro_id, collapse = ", ")
          interpro_description <- paste(interpro_description, collapse = "; ")
        } else if(length(interpro) == 0){
          interpro_id <- ""
          interpro_description <- ""
        }
      }
      c(description, ncbi_id, vectorbase_id, orthodb_id, uniprot_id, gos_id, 
        gos_description, interpro_id, interpro_description)
    }))
  colnames(annotation) <- c("description", "ncbi_id", "vectorbase_id", 
                            "orthodb_id", "uniprot_id", "GOs_id", 
                            "GOs_description", "interpro_id", 
                            "interpro_description")
  df <- cbind(df, annotation)
  return(df)
}

## Initialization =============================================================
# Initial variables
p1 <- 0.05 # First percentile
p2 <- 0.95 # Last percentile
gff_file <- "~/mosquito/data/refgenome/CulQui.gff"
abba.baba_table <- "~/mosquito/results/2019-02-23/abba_baba_slidingWindows.tsv"

# Load the data
gff <- read.table(gff_file, sep = "\t", comment.char = "#")
colnames(gff) <- c("chr", "source", "type", "start", "end", "score", "strand",
                   "cds", "atrribute")
gff <- gff[, -c(2, 6, 7, 8)]

f <- read.table(abba.baba_table)
f <- f[ ,-c(3, 5, 6, 7, 9, 10, 12)]

# percentile 5 and 95 in f1 - f2 (diff_f)
p <- list("05" = f[f$diff_f <= quantile(f$diff_f, p1), ],
          "95" = f[f$diff_f >= quantile(f$diff_f, p2), ])
## Some info
cat("\nThere are", as.character(nrow(p[["05"]])),
    "windows within the lowest 5% of the diff_f values")
cat("\nThere are", as.character(nrow(p[["95"]])),
    "windows within the highest 5% of the diff_f values")

## Annotation =================================================================
if(!dir.exists("annotation")){
  dir.create("annotation")
}
for(i in c("05", "95")){
  # Annotate
  ## Merge the annotation of the gff with the selected regions (windows)
  annot <- merge(p[[i]], gff, by = "chr")
  annot <- annot[annot$start > annot$pos &
                   annot$end < annot$pos + 1e5,]
  annot <- annot[annot$type == "gene", ]
  ## Retrieve more information on online databases
  annot <- annotate_culex(annot)
  
  # Save the results...
  ## ... annotated table
  write.table(annot, file = paste0("annotation/diff_f-percentile", i, ".tsv"), 
              sep = "\t", quote = F, row.names = F, col.names = T)
  ## ... lists of ncbi ids
  write.table(annot[, c("ncbi_id", "diff_f")], 
              file = paste0("annotation/ncbi_ids_p", i, ".tsv"), 
              sep = "\t", quote = F, row.names = F, col.names = F)
  ## ... lists of uniprot ids
  write.table(annot[annot$uniprot_id != "", c("uniprot_id", "diff_f")],
              file = paste0("annotation/uniprot_ids_p", i,".tsv"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  ## ... lists of GO ids
  write.table(unlist(strsplit(as.character(annot[annot$GOs_id != "", "GOs_id"]),
                              ", ")),
              file = paste0("annotation/GO_ids_p", i, ".tsv"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  cat("\nThere are", nrow(annot),
      "annotated genes within the windows in the percentile", i, 
      "of the diff_f values")
  cat("\nThe associated functions to them are these:\n")
  cat(paste0(annot$GOs_id[annot$GOs_id != ""], collapse = ", "))
  cat("\n\n")
}; rm(i, annot)

## Functional Analysis (1): Network Analysis -- STRING ========================
#' By using STRING, it is possible to know how the proteins in each list
#' interact. STRING returns the protein-protein interaction network and 
#' some statistics about the composition of the network
if(!dir.exists("string")){
  dir.create("string")
}
# Load the uniprot ids
for(i in c("05", "95")){
  uniprot <- read.table(paste0("annotation/uniprot_ids_p", i ,".tsv"),
                        header = F, sep = "\t",
                        col.names = c("id", "diff_f"), stringsAsFactors = F)
  uniprot <- paste0(uniprot$id, collapse = "%0d")
  
  # Retrieve the network image
  curl::curl_download(url = paste0(c("https://string-db.org/api/image/network?identifiers=",
                                     uniprot, "&species=7157"), collapse = ""),
                      destfile = paste0("string/string_network_p", i, ".png"),
                      mode = "w")
  
  # Retrieve the network as a table
  curl::curl_download(url = paste0(c("https://string-db.org/api/tsv/network?identifiers=",
                                     uniprot, "&species=7157"), collapse = ""),
                      destfile = paste0("string/string_network_p", i, ".tsv"),
                      mode = "w")
  
  # Retrieve statistics of the network structure
  curl::curl_download(url = paste0(c("https://string-db.org/api/tsv/ppi_enrichment?identifiers=",
                                     uniprot, "&species=7157"), collapse = ""),
                      destfile = paste0("string/string_network_p", i, "_stats.tsv"),
                      mode = "w")
}
uniprot_05 <- read.table("annotation/uniprot_ids_p05.tsv", header = F, sep = "\t",
                         col.names = c("id", "diff_f"), stringsAsFactors = F)
uniprot_05 <- paste0(uniprot_05$id, collapse = "%0d")
uniprot_95 <- read.table("annotation/uniprot_ids_p95.tsv", header = F, sep = "\t",
                         col.names = c("id", "diff_f"))
uniprot_95 <- paste0(uniprot_95$id, collapse = "%0d")

# Retrieve the network image
curl::curl_download(url = paste0(c("https://string-db.org/api/image/network?identifiers=",
                                   uniprot_05, "&species=7157"), collapse = ""),
                    destfile = "string/string_network_p05.png", mode = "w")
curl::curl_download(url = paste0(c("https://string-db.org/api/image/network?identifiers=",
                                   uniprot_95, "&species=7157"), collapse = ""),
                    destfile = "string/string_network_p95.png", mode = "w")

# Retrieve the network as a table
curl::curl_download(url = paste0(c("https://string-db.org/api/tsv/network?identifiers=",
                                   uniprot_05, "&species=7157"), collapse = ""),
                    destfile = "string/string_network_p05.tsv", mode = "w")
curl::curl_download(url = paste0(c("https://string-db.org/api/tsv/network?identifiers=",
                                   uniprot_95, "&species=7157"), collapse = ""),
                    destfile = "string/string_network_p95.tsv", mode = "w")

# Retrieve statistics of the network structure
curl::curl_download(url = paste0(c("https://string-db.org/api/tsv/ppi_enrichment?identifiers=",
                                   uniprot_05, "&species=7157"), collapse = ""),
                    destfile = "string/string_network_p05_stats.tsv", mode = "w")
curl::curl_download(url = paste0(c("https://string-db.org/api/tsv/ppi_enrichment?identifiers=",
                                   uniprot_95, "&species=7157"), collapse = ""),
                    destfile = "string/string_network_p95_stats.tsv", mode = "w")
rm(uniprot_05, uniprot_95)

## Functional Analysis (2): Over Representation Analysis (ORA) ================
# Creation of the DBi ORG package.
if(!dir.exists("annotation/CulQui")){
  library(AnnotationForge)
  dir.create("annotation/CulQui")
  makeOrgPackageFromNCBI(version = "0.1",
                         author = "Victor Moreno-Gonzalez <vicmogon@alumni.uv.es>",
                         maintainer = "Victor Moreno-Gonzalez <vicmogon@alumni.uv.es>",
                         outputDir = "./annotation/CulQui",
                         tax_id = "7176", # C. quinquefasciatus taxid
                         genus = "Culex",
                         species = "quinquefasciatus")
  install.packages("./annotation/CulQui/org.Cquinquefasciatus.eg.db", 
                   repos = NULL)
}
library(org.Cquinquefasciatus.eg.db)

# ORA
library(clusterProfiler)
if(!dir.exists("EnrichmentAnalysis")){
  dir.create("EnrichmentAnalysis")
}
for(i in c("05", "95")){
  # Load the data
  ncbi_ids <- read.table(paste0("annotation/ncbi_ids_p", i, ".tsv"),
                         header = F, sep = "\t",
                         col.names = c("id", "diff_f"), stringsAsFactors = F)
  ncbi_ids <- as.character(ncbi_ids$id)
  # Enrichment Analysis
  ego <- enrichGO(gene = ncbi_ids,
                  OrgDb = org.Cquinquefasciatus.eg.db,
                  keyType = "ENTREZID",
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  # Plot results
  if(nrow(ego@result) > 0){
    png(filename = paste0("EnrichmentAnalysis/diff-f_cnetplot_", i, ".png"),
        width = 30, height = 30, units = "cm", res = 333)
    cnetplot(ego, font.size = 5)
    dev.off()
  }
  # Save the results
  write.table(ego@result, file = paste0("EnrichmentAnalysis/diff-f_p", i, ".tsv"),
              quote = F, sep = "\t", col.names = T)
}; rm(i, ncbi_ids, ego)

