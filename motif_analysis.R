###setup ####
libs <- c("dplyr",
          "pheatmap",
          "tibble",
          "readxl",
          "openxlsx",
          "tidyr",
          "ggplot2",
          "biomaRt",
          'magrittr',
          'ggseqlogo')
lapply(libs, require, character.only=T)

### input gene list in official gene symbol format ####
# input is the list of 1032 differentially expressed genes
input <- readRDS("input.rds")
genelist <- list(input$Gene)

### pull in fpkm values for all genes ####
fpkm <- readRDS("fpkm.rds")

### extract principle isoform based on APPRIS ####
mart <- useMart("ensembl")
listDatasets(mart)
mart <- useDataset("mmusculus_gene_ensembl", mart)
att <- listAttributes(mart)
f <- listFilters(mart)
atributes <- c("external_gene_name",'ensembl_gene_id','transcript_appris',
               'ensembl_transcript_id','refseq_mrna','transcription_start_site',
               'chromosome_name', 'strand', "start_position", "end_position")
filters <- "ensembl_gene_id"
values <- genelist
output <- getBM(attributes=atributes, filters=filters, values=values, mart=mart)
ens2gene <- getBM(attributes=c("external_gene_name","ensembl_gene_id"), mart=mart) %>% unique()
output_principle <- output[grepl('principal', output$transcript_appris),]
output_principle$strand[output_principle$strand<0] <- "-"
output_principle$strand[output_principle$strand>0] <- "+"
output_principle <- output_principle %>% left_join(input[,c(1,4)], by = c("ensembl_gene_id"="Gene"))
colnames(output_principle)[11] <- "score"
output_principle <- output_principle[grepl("NM", output_principle$refseq_mrna),] %>% 
  unique() %>%
  dplyr::select(-c(refseq_mrna)) %>%
  unique()


### get pos and neg strands to calculate ROIs ####
output_principle_pos <- output_principle %>%
  filter(strand == '+') %>%
  mutate("tenkbup"=transcription_start_site-10000) %>%
  mutate("fivekbup"=transcription_start_site-5000) %>%
  mutate("onekbup"=transcription_start_site-1000)%>%
  mutate("thdown"=transcription_start_site+200)
output_principle_pos <- output_principle_pos[order(output_principle_pos$external_gene_name, 
                                                   output_principle_pos$transcription_start_site, 
												   decreasing=F),]

output_principle_pos <- output_principle_pos[!duplicated(output_principle_pos$external_gene_name),]

output_principle_neg <- output_principle %>%
  filter(strand == '-') %>%
  mutate("tenkbup"=transcription_start_site+10000) %>%
  mutate("fivekbup"=transcription_start_site+5000) %>%
  mutate("onekbup"=transcription_start_site+1000) %>%
  mutate("thdown"=transcription_start_site-200)
output_principle_neg <- output_principle_neg[order(output_principle_neg$external_gene_name, 
                                                   output_principle_neg$transcription_start_site, 
												   decreasing=T),]

output_principle_neg <- output_principle_neg[!duplicated(output_principle_neg$external_gene_name),]

### figure out if we are missing any genes from the input list
final_genes <- c(output_principle_pos$ensembl_gene_id, output_principle_neg$ensembl_gene_id)
start_genes <- input$Gene
missing_genes <- setdiff(start_genes, final_genes)

missing <- output[grep(paste(missing_genes, collapse="|"),output$ensembl_gene_id),]
missing$strand[missing$strand<0] <- "-"
missing$strand[missing$strand>0] <- "+"
missing <- missing %>% left_join(input[,c(1,4)], by = c("ensembl_gene_id"="Gene"))
colnames(missing)[11] <- "score"
missing <- missing[order(missing$ensembl_transcript_id, decreasing=F),]
missing <- missing[!duplicated(missing$external_gene_name),]

missing_pos <- missing %>%
  filter(strand == '+') %>%
  mutate("onekbup"=start_position-1000) %>%
  mutate("thdown"=start_position+200)
missing_neg <- missing %>%
  filter(strand == '-') %>%
  mutate("onekbup"=end_position+1000) %>%
  mutate("thdown"=end_position-200)

### make stranded bed files ####


onekbupthdown_bed_pos_mis <- missing_pos %>% 
  as.data.frame() %>% 
  dplyr::select(chromosome_name, onekbup, thdown, ensembl_gene_id, score, strand) %>% 
  set_colnames(c("chrom","start","stop","name","score","strand"))


onekbupthdown_bed_neg_mis <- missing_neg %>% 
  as.data.frame() %>% 
  dplyr::select(chromosome_name, thdown, onekbup, ensembl_gene_id, score, strand) %>% 
  set_colnames(c("chrom","start","stop","name","score","strand"))


onekbupthdown_bed_pos <- output_principle_pos %>% 
  as.data.frame() %>% 
  dplyr::select(chromosome_name, onekbup, thdown, ensembl_gene_id, score, strand) %>% 
  set_colnames(c("chrom","start","stop","name","score","strand"))

onekbupthdown_bed_neg <- output_principle_neg %>% 
  as.data.frame() %>% 
  dplyr::select(chromosome_name, thdown, onekbup, ensembl_gene_id, score, strand) %>% 
  set_colnames(c("chrom","start","stop","name","score","strand"))

### combine standed bed files ####


onekbupthdown_bed <- rbind(onekbupthdown_bed_pos,onekbupthdown_bed_neg, onekbupthdown_bed_pos_mis, onekbupthdown_bed_neg_mis) %>% mutate(ROI=paste0(start,"-",stop,":",strand))


write.table(select(onekbupthdown_bed, c(chrom, start, stop, name, score, strand)), file="1kbup200down.bed", sep = "\t", quote = F, col.names = F, row.names = F)

########################
# Now that BED files are created, switch to bash.

# 1) convert bed file to fasta
```
bedtools getfasta -fi Mus_musculus.GRCm38.dna.primary_assembly.fa -bed 1kbup200down.bed > 1KB_up_TSS.fa
```
# 2) run Ame on the fasta
```
ame --o 1KB_up_ame_HOC --control shuffle 1KB_up_TSS.fa HOCOMOCOv11_full.meme
```
# 3) parse Ame output below

#######################



### parse AME sequences.tsv ####
#Sox2 = M6129_1.02
#Foxa2 = M1879_1.02
#Foxa1 = M1965_1.02
temp <- read.table("1KB_up_ame_HOC\sequences.tsv", sep="\t", stringsAsFactors = F, header =T) %>%
  separate(motif_ID, c("motif"), sep="_")
temp_sub <- temp[grep("SOX2|FOXA1|FOXA2", temp$motif),]
roi2_gene <- input[,c(1,2)] %>% unique()
temp_sub_annot <- temp_sub %>% left_join(roi2_gene, by = c("seq_ID" = "Gene")) %>%
  na.omit() %>% unique()
l <- list("Sox2" = temp_sub_annot[temp_sub_annot$motif == "SOX2",],
          "Foxa2" = temp_sub_annot[temp_sub_annot$motif == "FOXA2",],
          "Foxa1" = temp_sub_annot[temp_sub_annot$motif == "FOXA1",])
openxlsx::write.xlsx(l, "1kbup200down_ame_hoco_DIM_FoxSox_hits.xlsx")

ame <- read.table("results/1kbup200down_ame_hoco_dim/ame.tsv",sep="\t",stringsAsFactors = F, header = T) %>%
  separate(motif_ID, c("motif"),sep="_")
ame$motif <- lapply(tolower(ame$motif), simpleCap) %>% as.character()
ame <- ame %>% left_join(ens2gene, by=c("motif"="external_gene_name")) %>%
  left_join(fpkm, by=c("ensembl_gene_id"="Gene")) %>%
  dplyr::select(c(-1,-2,-4, -22)) %>% dplyr::select(c(1,15:18,2:14)) %>%
  mutate("ratio" = TP/FP) %>% dplyr::select(c(1:5,19,6:18))
l <- list("ame" = ame,
          "filtered"= ame[ame$dim > 5 | ame$bright > 5,])
openxlsx::write.xlsx(l, "1kbup200down_ame_DIM_motif_hits.xlsx")

simpleCap <- function(x) {
  tmp <- strsplit(x, " ")[[1]]
  paste(toupper(substring(tmp, 1, 1)), 
		substring(tmp, 2),
        collapse = " ", 
		sep = "")
}