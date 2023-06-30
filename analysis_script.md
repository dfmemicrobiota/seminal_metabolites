# Impact of semen microbiota on the composition of seminal plasma![image](https://github.com/dfmemicrobiota/seminal_metabolites/assets/75491874/03d6dc6e-b40a-40d7-8cde-e548915834a4)

### Description
Male infertility is a significant concern contributing to the growing utilization of assisted reproductive technologies (ART). Environmental and lifestyle factors have been implicated, and the impact of the male genital microbiota on fertility is an unexplored area of research. This study aimed to investigate the association between seminal microbiota and sperm physiology through the integration of microbiota profiling and metabolomics.

## Import libraries

```
library(htmltools)
library(utils)
library(BiocManager)
library(DESeq2)
library(ggplot2)
library(lefser)
library(EnhancedVolcano)
library(dplyr)
library(corrplot)
library(tidyverse)
library(pheatmap)
library(lefser)
library(reshape)
library(plyr)
library(RColorBrewer)
library(pals)
library(ggpubr)
library(hmdbQuery)
library(XML)
library(matrixStats)
library(gridExtra)
library(grid)
library(scales)
```

## Import microbiota data (phyloseq object)
```
otu_table<-otu_table(seqtab.nochim, taxa_are_rows=F)
tax_table_taxa<-tax_table(taxa)

write.table(otu_table, "~/sperm_metabolites/otu_table.csv", sep="\t")
write.table(samdf, "~/sperm_metabolites/samdf.csv", sep="\t")
write.table(tax_table_taxa, "~/sperm_metabolites/tax_table_taxa.csv", sep="\t")

setwd("~/Google Drive/PROJECTS/SIC metabolomics/input")

otu_table<-read.csv("~/sperm_metabolites/otu_table.csv", sep="\t")
samdf<-read.csv("~/sperm_metabolites/samdf.csv", sep="\t")
tax_table_taxa<-read.csv("~/sperm_metabolites/tax_table_taxa.csv", sep="\t")

```

## Import metabolite data and subset

```
setwd("~/sperm_metabolites/")

#import metabolite data
countData_adjusted <- read.csv("~/sperm_metabolites/input_volnorm_adjusted_updated.csv", header = TRUE, sep = ",")

```

## Parse all HMDB xml file (Python)

### Script used to parse the xml file
```
from io import StringIO
from lxml import etree
import csv
def hmdbextract(name, file):
  ns = {'hmdb': 'http://www.hmdb.ca'}
  context = etree.iterparse(name, tag='{http://www.hmdb.ca}metabolite')
  csvfile = open(file, 'w')
  fieldnames = ['accession', 'synonym1', 'synonym2', 'synonym3', 'synonym4', 'synonym5', 'synonym6','iupac_name', 'name', 'chemical_formula', 'InChIKey', 'cas_registry_number', 'smiles', 'drugbank','chebi_id', 'pubchem', 'phenol_explorer_compound_id','food','knapsack', 'chemspider', 'kegg', 'meta_cyc','bigg','metlin_id','pdb_id', 'logpexp','kingdom',  'direct_parent', 'super_class', 'class', 'sub_class', 'molecular_framework']
  writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
  writer.writeheader()
  for event, elem in context:

    accession = elem.xpath('hmdb:accession/text()', namespaces=ns)[0]
    
    try:
        synonym1 = elem.xpath('hmdb:synonyms/hmdb:synonym/text()', namespaces=ns)[0]
    except:
        synonym1 = 'NA'
        
    try:
        synonym2 = elem.xpath('hmdb:synonyms/hmdb:synonym/text()', namespaces=ns)[1]
    except:
        synonym2 = 'NA'
        
    try:
        synonym3 = elem.xpath('hmdb:synonyms/hmdb:synonym/text()', namespaces=ns)[2]
    except:
        synonym3 = 'NA'
        
    try:
        synonym4 = elem.xpath('hmdb:synonyms/hmdb:synonym/text()', namespaces=ns)[3]
    except:
        synonym4 = 'NA'
    
    try:
        synonym5 = elem.xpath('hmdb:synonyms/hmdb:synonym/text()', namespaces=ns)[4]
    except:
        synonym5 = 'NA'
    
    try:
        synonym6 = elem.xpath('hmdb:synonyms/hmdb:synonym/text()', namespaces=ns)[5]
    except:
        synonym6 = 'NA'

    try:
        iupac_name = elem.xpath('hmdb:iupac_name/text()', namespaces=ns)[0].encode('utf-8')
    except:
        iupac_name = 'NA'
    name = elem.xpath('hmdb:name/text()', namespaces=ns)[0].encode('utf-8')
    try:
        chemical_formula = elem.xpath('hmdb:chemical_formula/text()', namespaces=ns)[0]
    except:
        chemical_formula = 'NA'
    try:
        inchikey = elem.xpath('hmdb:inchikey/text()', namespaces=ns)[0]
    except:
        inchikey = 'NA'
    try:
        cas_registry_number = elem.xpath('hmdb:cas_registry_number/text()', namespaces=ns)[0]
    except:
        cas_registry_number = 'NA'
    try:
        smiles = elem.xpath('hmdb:smiles/text()', namespaces=ns)[0]
    except:
        smiles = 'NA'
    try:
        drugbank = elem.xpath('hmdb:drugbank_id/text()', namespaces=ns)[0]
    except:
        drugbank = 'NA'
    try:
        chebi_id = elem.xpath('hmdb:chebi_id/text()', namespaces=ns)[0]
    except:
        chebi_id = 'NA'
    try:
        pubchem = elem.xpath('hmdb:pubchem_compound_id/text()', namespaces=ns)[0]
    except:
        pubchem = 'NA'
    try:
        phenol_explorer_compound_idt = elem.xpath('hmdb:phenol_explorer_compound_id/text()', namespaces=ns)[0]
    except:
        phenol_explorer_compound_id = 'NA'
    try:
        food = elem.xpath('hmdb:foodb_id/text()', namespaces=ns)[0]
    except:
        food = 'NA'
    try:
        knapsack = elem.xpath('hmdb:knapsack_id/text()', namespaces=ns)[0]
    except:
        knapsack = 'NA'
    try:
        chemspider = elem.xpath('hmdb:chemspider_id/text()', namespaces=ns)[0]
    except:
        chemspider = 'NA'
    try:
        kegg = elem.xpath('hmdb:kegg_id/text()', namespaces=ns)[0]
    except:
        kegg = 'NA'
    try:
        meta_cyc = elem.xpath('hmdb:meta_cyc_id/text()', namespaces=ns)[0]
    except:
        meta_cyc = 'NA'
    try:
        bigg = elem.xpath('hmdb:bigg_id/text()', namespaces=ns)[0]
    except:
        bigg = 'NA'
    try:
        metlin_id = elem.xpath('hmdb:metlin_id/text()', namespaces=ns)[0]
    except:
        metlin_id = 'NA'
    try:
        pdb_id = elem.xpath('hmdb:pdb_id/text()', namespaces=ns)[0]
    except:
        pdb_id = 'NA'
    try:
        logpexp = elem.xpath('hmdb:experimental_properties/hmdb:property[hmdb:kind = "logp"]/hmdb:value/text()', namespaces=ns)[0]
    except:
        logpexp = 'NA'
    try:
        kingdom = elem.xpath('hmdb:taxonomy/hmdb:kingdom/text()', namespaces=ns)[0]
    except:
        kingdom = 'NA'
    try:
        direct_parent = elem.xpath('hmdb:taxonomy/hmdb:direct_parent/text()', namespaces=ns)[0]
    except:
        direct_parent = 'NA'
    try:
        super_class = elem.xpath('hmdb:taxonomy/hmdb:super_class/text()', namespaces=ns)[0]
    except:
        super_class = 'NA'
    try:
        classorg = elem.xpath('hmdb:taxonomy/hmdb:class/text()', namespaces=ns)[0]
    except:
        classorg = 'NA'
    try:
        sub_class = elem.xpath('hmdb:taxonomy/hmdb:sub_class/text()', namespaces=ns)[0]
    except:
        sub_class = 'NA'
    try:
        molecular_framework = elem.xpath('hmdb:taxonomy/hmdb:molecular_framework/text()', namespaces=ns)[0]
    except:
        molecular_framework = 'NA'

    writer.writerow({'accession': accession, 'synonym1': synonym1,'synonym2': synonym2,'synonym3': synonym3,'synonym4': synonym4,'synonym5': synonym5,'synonym6': synonym6, 'iupac_name': iupac_name, 'name': name, 'chemical_formula': chemical_formula, 'InChIKey': inchikey, 'cas_registry_number': cas_registry_number, 'smiles': smiles,'drugbank': drugbank,'chebi_id': chebi_id,'pubchem': pubchem,'phenol_explorer_compound_id':phenol_explorer_compound_id, 'food': food,'knapsack': knapsack, 'chemspider': chemspider,'kegg': kegg, 'meta_cyc': meta_cyc, 'bigg':bigg, 'metlin_id': metlin_id, 'pdb_id':pdb_id,'logpexp':logpexp, 'kingdom': kingdom, 'direct_parent': direct_parent, 'super_class': super_class, 'class': classorg, 'sub_class': sub_class, 'molecular_framework': molecular_framework})
    # It's safe to call clear() here because no descendants will be
    # accessed
    elem.clear()
# Also eliminate now-empty references from the root node to elem
    for ancestor in elem.xpath('ancestor-or-self::*'):
        while ancestor.getprevious() is not None:
            del ancestor.getparent()[0]
  del context
  return;
```


### Run the script used to parse the xml file
```

hmdbextract('HMDBmissing.xml','missing.csv')

#remove all "<?xml version="1.0" encoding="UTF-8"?>" except the first - add "<hmdb xmlns="http://www.hmdb.ca">" at the second line - add "</hmdb>" at the very end
```

### Process the data and classify metabolites according to HMDB data

```
all_hmdb<- read.csv('all_hmdb.csv', header = TRUE, sep = ",")

colnames(all_hmdb)

all_hmdb$synonyms<-paste(all_hmdb$synonym1,all_hmdb$synonym2,all_hmdb$synonym3,all_hmdb$synonym4,all_hmdb$synonym5,all_hmdb$synonym6,sep=" ")

write.table(all_hmdb, "~/sperm_metabolites/all_hmdb.csv", sep=",")

present<-data.frame(matrix(NA, nrow = 0, ncol = 8))

#Find matches between metab of the countData_adjusted and metabolite names in the database
for (i in 1:length(countData_adjusted$metab)) {

tryCatch({ 
  for (z in 1:length(all_hmdb$name)) {
    
    if(grepl(countData_adjusted$metab[i], all_hmdb$name[z]) == TRUE){
    
      present<-rbind(present,c(countData_adjusted$metab[i],all_hmdb$accession[z],all_hmdb$name[z],
                               all_hmdb$kingdom[z],all_hmdb$direct_parent[z],all_hmdb$super_class[z],
                               all_hmdb$class[z],all_hmdb$sub_class[z],all_hmdb$molecular_framework[z]))
      }
    }
  
  }, error = function(e) -999)
}

colnames(present)<-c("metab","name","kingdom","direct_parent","super_class","class","sub_class","molecular_framework")
present<-unique(present)

#Repeat the search for entries present in synonyms
for (i in 1:length(countData_adjusted$metab)) {
  
  tryCatch({ 
    for (z in 1:length(all_hmdb$synonyms)) {
      
      if(grepl(countData_adjusted$metab[i], all_hmdb$synonyms[z]) == TRUE){
        
        present<-rbind(present,c(countData_adjusted$metab[i],all_hmdb$accession[z],all_hmdb$name[z],
                                 all_hmdb$kingdom[z],all_hmdb$direct_parent[z],all_hmdb$super_class[z],
                                 all_hmdb$class[z],all_hmdb$sub_class[z],all_hmdb$molecular_framework[z]))
      }
    }
  }, error = function(e) -999)
}

#rename the table columns
colnames(present)<-c("metab","accession","name","kingdom","direct_parent","super_class","class","sub_class","molecular_framework")
present<-unique(present)

#Save output and manually curate file
write.table(present, "~/sperm_metabolites/present.csv", sep=",")
```

## Load the data used for the analysis

```
setwd("~/sperm_metabolites/")

countData_adjusted <- read.csv("~/sperm_metabolites/input_volnorm_adjusted_updated.csv", header = TRUE, sep = ",")

present_curated<-read.csv("~/sperm_metabolites/present_curated.csv", sep=",",header=T)
present<-unique(present_curated)
colnames(present_curated)<-c("metab","accession","name","kingdom","direct_parent","super_class","class","sub_class","molecular_framework")

metab_classes<-left_join(present_curated,countData_adjusted,by="metab")
write.table(metab_classes, "~/sperm_metabolites/metab_classes.csv", sep=",")

#Import final file, manually curated
metab_classes<-read.csv("~/sperm_metabolites/metab_classes.csv", sep=",",header=T)

#Create files for each metabolite classification group
metab_count_data<-metab_classes[,10:28]
metab_count<-metab_count_data
metab_count$taxonomy<-metab_classes$metab

metab_count_direct_parent<-metab_count_data
metab_count_direct_parent$taxonomy<-metab_classes$direct_parent
metab_count_direct_parent<-ddply(metab_count_direct_parent, "taxonomy", numcolwise(sum))

metab_count_super_class<-metab_count_data
metab_count_super_class$taxonomy<-metab_classes$super_class
metab_count_super_class<-ddply(metab_count_super_class, "taxonomy", numcolwise(sum))

metab_count_class<-metab_count_data
metab_count_class$taxonomy<-metab_classes$class
metab_count_class<-ddply(metab_count_class, "taxonomy", numcolwise(sum))

metab_count_sub_class<-metab_count_data
metab_count_sub_class$taxonomy<-metab_classes$sub_class
metab_count_sub_class<-ddply(metab_count_sub_class, "taxonomy", numcolwise(sum))
```

## Metabolite description - Figure 1


### Show only top10 
```
#SWITCH
input<-metab_count_sub_class
input[is.na(input)] = "NA"

input$metab<-NULL
input<-ddply(input,"taxonomy",numcolwise(sum))
length.db<-dim(input)[1]

#absolute values
row.names(input)<-input$taxonomy
input$taxonomy<-NULL
input<-rowSums(input)
input<-rev(sort(input))

input_<-c(input[1:4],input[6:11],"Others"=(sum(input[12:length.db])+input["NA"]))
names(input_)[11]<-"Others"

input_<-melt(input_)
input_<-data.frame(input_,"taxonomy"=rownames(input_))

#Define graph colors
# library(randomcoloR)
# n <- dim(input_)[1]
# mycolors <- distinctColorPalette(n)

mycolors<-c("#A8E750","#DD60B6", "#AF4EE5", "#76D9D4", "#9B86DF", "#E0ACCB", "#84DD95", "#DC7B62", "#D5D9C4", "#89AED6","#D6C970")

pie_chart<-ggplot(input_, aes(x="",y=value,fill=taxonomy)) +
              geom_bar(stat="identity", width=1, color="white") +
              coord_polar("y", start=0) +
              scale_fill_manual(values = mycolors) +
              theme(legend.text=element_text(size=rel(2.5))) +
              theme_void() +
              guides(fill=guide_legend(title="Metabolite subclass"))


# Create barplot for the most abundant metabolites

# SWITCH
input<-metab_count
#input$metab<-rownames(input)

# Convert data to melted format
table_melt <- melt(input, by="taxonomy")
table_melt$variable<-NULL

# Calculate mean and standard deviation for each metabolite
summary_table<-ddply(table_melt, c("taxonomy"), summarise,
                     mean = mean(value), sd = sd(value),
                     sem = sd(value)/sqrt(length(value)))

# Order categories by mean value
summary_table <- as.data.frame(summary_table[rev(order(summary_table$mean)),])

#SWITCH
summary_table_top10<-summary_table[1:10,]
summary_table_top20<-summary_table[1:20,]
summary_table_top50<-summary_table[1:50,]

input<-summary_table_top20

color<-as.data.frame(cbind(metab_classes$metab,metab_classes$sub_class))
colnames(color)<-c("taxonomy","sub_class")
color<-left_join(input,color,by="taxonomy")
color[is.na(color)] ="Others"

colors_piechart<-data.frame(cbind(as.data.frame(ggplot_build(pie_chart)[[1]])$fill,sort(pie_chart$data$taxonomy)))
colnames(colors_piechart)<-c("code","sub_class")

colors_piechart<-as.data.frame(ggplot_build(pie_chart)[[1]])
colors_piechart<-colors_piechart[order(colors_piechart$group),]
colors_piechart<-data.frame(cbind(colors_piechart$fill,sort(pie_chart$data$taxonomy)))
colnames(colors_piechart)<-c("code","sub_class")

colors<- vector()

  for(i in 1:nrow(color)){
    for(z in 1:nrow(colors_piechart)){
    
    if(color$sub_class[i] == colors_piechart$sub_class[z]){
      
      colors<-c(colors,colors_piechart$code[z])
    }
   }
  }


# Create barplot with error bars
bar_plot <- ggplot(summary_table_top20, aes(x=reorder(taxonomy,-mean),y=mean,fill=taxonomy)) +
  geom_bar(stat = "identity",fill = colors) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  labs(x = "Metabolite", y = "Absolute quantification") +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  theme_minimal() +
  coord_cartesian(ylim=c(100000,100000000)) +
  theme_classic() +
  theme(text = element_text(size=12),axis.text.x = element_text(angle = 45, hjust=1))

# Print barplot to screen
print(bar_plot)
bar_plot
```

## Volcano plots - Figure 2

```
#switch classes
input<-metab_count

input<-ddply(input,"taxonomy",numcolwise(sum))
length.db<-dim(input)[1]
if(is.na(input$taxonomy[length.db])) {
 input<-input[-(length.db),]
 }  
colnames(input)[1]<-"metab"

#all_metabolites
input<-metab_count
rownames(input)<-input$taxonomy
input<-input %>% relocate(taxonomy)

#round the data to avoid decimals
input[,-1]<-round(input[,-1])
head(input)

#import metadata
metadata <- read.csv("~/sperm_metabolites/input_metadata.csv", header = TRUE, sep = ",")
#head(metadata)

clinicalData <- read.csv("~/sperm_metabolites/clinical_data.csv", header = TRUE, sep = ",")
#head(clinicalData)

#create a DESEQ2 object
dds_microbiota <- DESeqDataSetFromMatrix(countData=input,colData=metadata, design=~microbiota, tidy = TRUE)

#run DESEQ2 analysis
dds_microbiota <- DESeq(dds_microbiota)
res_microbiota <- results(dds_microbiota)
head(res_microbiota, tidy=TRUE) #let's look at the results table

# library(EnhancedVolcano)
# library(magrittr)
EnhancedVolcano(res_microbiota, lab = rownames(res_microbiota),
                x = 'log2FoldChange', y = 'pvalue',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,FCcutoff = 2.0,
                pointSize = 2.0,labSize = 3.5,colAlpha = 1,
                legendPosition = 'right',legendLabSize = 12,legendIconSize = 4.0,
                drawConnectors = TRUE,widthConnectors = 0.75,
                axisLabSize = 12,title = "",subtitle="")

#spermiogram

#clinical variables
dds_spermiogram <- DESeqDataSetFromMatrix(countData=input, colData=metadata,design=~spermiogram, tidy = TRUE)
dds_spermiogram <- DESeq(dds_spermiogram)
res_spermiogram<-results(dds_spermiogram)

#concentration
dds_conc <- DESeqDataSetFromMatrix(countData=input,colData=metadata, design=~conc, tidy = TRUE)
dds_conc <- DESeq(dds_conc)
res_conc<-results(dds_conc)

#count
dds_count <- DESeqDataSetFromMatrix(countData=input,colData=metadata,design=~count, tidy = TRUE)
dds_count <- DESeq(dds_count)
res_count<-results(dds_count)

#progressive motility
dds_prog_mot <- DESeqDataSetFromMatrix(countData=input,colData=metadata,design=~prog_mot, tidy = TRUE)
dds_prog_mot <- DESeq(dds_prog_mot)
res_prog_mot<-results(dds_prog_mot)

#tot motility
dds_tot_mot <- DESeqDataSetFromMatrix(countData=input,colData=metadata,design=~tot_mot, tidy = TRUE)
dds_tot_mot <- DESeq(dds_tot_mot)
res_tot_mot<-results(dds_tot_mot)

#morphology
dds_morph <- DESeqDataSetFromMatrix(countData=input,colData=metadata,design=~morph, tidy = TRUE)
dds_morph <- DESeq(dds_morph)
res_morph<-results(dds_morph)

#all
dds_all <- DESeqDataSetFromMatrix(countData=input,colData=metadata,design=~prog_mot+tot_mot, tidy = TRUE)
dds_all <- DESeq(dds_all)
res_all<-results(dds_all)

###CREATE THE VOLCANO PLOT - enhanced volcano###

input_plot<-res_spermiogram
input_plot<-res_prog_mot
input_plot<-res_tot_mot
input_plot<-res_conc
input_plot<-res_count
input_plot<-res_morph
input_plot<-res_all
# library(EnhancedVolcano)
# library(magrittr)

plot <- EnhancedVolcano(input_plot, lab = rownames(input_plot),
            x = 'log2FoldChange', y = 'pvalue',
            xlab = bquote(~Log[2]~ 'fold change'),
            pCutoff = 0.05,FCcutoff = 2.0,
            pointSize = 2.0,labSize = 5,colAlpha = 1,
            legendPosition = 'right',legendLabSize = 12,legendIconSize = 4.0,
            drawConnectors = TRUE,widthConnectors = 0.75,
            axisLabSize = 12,
            title = unlist(strsplit(input_plot@elementMetadata$description[3], ": "))[2],
            subtitle="")


# save plot if needed
pdf("volcanoplot.pdf", width = 14,height = 6)

EnhancedVolcano(input_plot, lab = rownames(input_plot),
                x = 'log2FoldChange', y = 'pvalue',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,FCcutoff = 2.0,
                pointSize = 4.0,labSize = 3.5,colAlpha = 1,
                legendPosition = 'right',legendLabSize = 10,legendIconSize = 4.0,
                drawConnectors = TRUE,widthConnectors = 0.75)
dev.off()
```

## Statistical analysis spermiogram parameters vs metabolites - Figure 3

```
cbind.fill <- function(...){
    nm <- list(...) 
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

subsetList <- function(myList, elementNames) {
    lapply(elementNames, FUN=function(x) myList[[x]])
}

input<-metab_count

input<-ddply(input,"taxonomy",numcolwise(sum))
length.db<-dim(input)[1]
input$taxonomy[length.db]<-"NA"
rownames(input)<-input$taxonomy
input$taxonomy<-NULL
input<-t(input)
View(input)

metadata <- read.csv("~/sperm_metabolites/input_metadata.csv")
metadata<-subset(metadata, select = c(id,microbiota, spermiogram, conc, count,tot_mot,prog_mot,morph))
rownames(metadata)<-metadata$id
metadata$id<-NULL
view(metadata)

#perform w-test
results_wtest<-data.frame(matrix(NA, nrow = 0, ncol = 3))

for (z in 1:ncol(metadata)) {

  for (i in 1:ncol(input)) {
  
  test<-as.data.frame(cbind(metadata[,z],input[,i]))

  data_table<-data.frame(cbind.fill(V1=as.numeric(unlist(test[test$V1 == unique(metadata[,z])[1],][2])),
                            V2=as.numeric(unlist(test[test$V1 == unique(metadata[,z])[2],][2]))))
  
  res_wtest<-wilcox.test(data_table$X1,data_table$X2)
  
  results_wtest<-rbind(results_wtest,c(colnames(metadata)[z],colnames(input)[i],res_wtest$p.value))
  colnames(results_wtest)<-c("metadata","metab","p-value")
  }
}

results_wtest.<-subset(results_wtest,results_wtest$`p-value`<0.05)
View(results_wtest.)

#Select only acyl carnitines for spermiogram analysis

acylcarnitines<-subset(metab_classes,metab_classes$direct_parent=="Acyl carnitines")
acylcarnitines<-acylcarnitines[,1:3]

results_wtest.<-left_join(acylcarnitines,results_wtest.,by="metab")

plots<-list()

for (z in 1:nrow(results_wtest.)) {

# results_wtest.$metadata[z]
# results_wtest.$metab[z]

test<-as.data.frame(cbind(metadata[,results_wtest.$metadata[z]],as.numeric(unlist(input[,results_wtest.$metab[z]]))))
test$V2<- sapply(test$V2, as.numeric)
#view(test)

plots[[z]]<-ggplot(test, aes(x=V1, y=V2)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  ylab(results_wtest.$metab[z]) +
  xlab("") +
  ggtitle(paste0(results_wtest.$metadata[z]," ~ ",results_wtest.$metab[z])) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)),
                breaks = c(100, 1000, 10000, 100000, 1000000)) +
  theme_bw() 
  
}

plots_spermiogram<-subsetList(plots, which(results_wtest.$metadata=="spermiogram"))
n <- length(plots_spermiogram)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(plots_spermiogram, ncol=4))

#Plots
my_comparisons <- list(c("abnormal", "normal"))

plots_spermiogram[[1]]+ 
  stat_compare_means(comparisons = my_comparisons,label = "p.format", method = "wilcox.test",label.y = log10(max(plots_spermiogram[[2]]$data$V2))) + 
  scale_y_log10(breaks = c(10000, 100000, 1000000),limits =c(10000,1000000),labels = trans_format("log10", math_format(10^.x))) + theme_minimal()
              
plots_spermiogram[[2]]+ 
  stat_compare_means(comparisons = my_comparisons,label = "p.format", method = "wilcox.test",label.y = log10(max(plots_spermiogram[[3]]$data$V2))) + 
  scale_y_log10(breaks = c(10000, 100000, 1000000),limits =c(10000,1000000),labels = trans_format("log10", math_format(10^.x))) + theme_minimal()

plots_spermiogram[[3]]+ 
  stat_compare_means(comparisons = my_comparisons,label = "p.format", method = "wilcox.test",label.y = log10(max(plots_spermiogram[[4]]$data$V2))) + 
  scale_y_log10(breaks = c(10000, 100000, 1000000),limits =c(10000,1000000),labels = trans_format("log10", math_format(10^.x))) + theme_minimal()

plots_spermiogram[[4]]+ 
  stat_compare_means(comparisons = my_comparisons,label = "p.format", method = "wilcox.test",label.y = log10(max(plots_spermiogram[[5]]$data$V2))) + 
  scale_y_log10(breaks = c(10000, 100000, 1000000),limits =c(10000,1000000),labels = trans_format("log10", math_format(10^.x))) + theme_minimal()
```

## Principal component analysis - Figure 4

```
#PCE for metabolites

#install.packages(c("FactoMineR", "factoextra"))

library("FactoMineR")
library("factoextra")

#If subclasses
metabolites<-metab_count
metabolites<-ddply(metabolites,"taxonomy",numcolwise(sum))
length.db<-dim(metabolites)[1]
if(is.na(metabolites$taxonomy[length.db])) {
  metabolites<-metabolites[-(length.db),]
}  

#For all metabolites
metabolites<-metab_count
row.names(metabolites)<-metabolites$taxonomy
metabolites$taxonomy<-NULL
metabolites<-t(metabolites)

metadata <- read.csv("~/sperm_metabolites/input_metadata.csv")
rownames(metadata)<-metadata$id
metadata$id<-NULL


s.metabolites<-scale(metabolites)
res.pca <- prcomp(s.metabolites)

#Microbiota information
microbiota <- read.csv("~/sperm_metabolites/microbiota_data.csv")
rownames(microbiota)<-microbiota$id
microbiota$id<-NULL
microbiota<-t(microbiota)

prevotella<-microbiota[,1]
lactobacillus<-microbiota[,3]

PvsL<-data.frame("prevotella"=prevotella,"lactobacillus"=lactobacillus)

PvsL$abundance<-PvsL$prevotella-PvsL$lactobacillus

for(n in 1:length(PvsL$abundance)){
  if(PvsL$abundance[n]<0){
    PvsL$microbiota_type[n]="Lactobacillus-type"
  }else{
    PvsL$microbiota_type[n]="Prevotella-type"
  }
}

PvsL$abundance<-abs(PvsL$abundance)


pca_plot<-fviz_pca_ind(res.pca, pointsize = PvsL$abundance,
             pointshape = 21, fill = PvsL$microbiota_type,
             repel = TRUE, # Avoid text overlapping (slow if many points)
              title= "", mean.point=F, label = "none",
             addEllipses = TRUE, ellipse.level = 0.95) + theme_classic()

pdf("~/sperm_metabolites/output/pca_microbiota.pdf", width =6,height = 4,bg = "transparent")
pca_plot
dev.off()




#spermiogram
pca_plot_spermiogram<-fviz_pca_ind(res.pca, pointsize = 3,
                       pointshape = 21, fill = metadata$spermiogram,
                       repel = TRUE, # Avoid text overlapping (slow if many points)
                       label = "none", title= "normal vs abnormal spermiogram",
                       mean.point=F,
                       addEllipses = TRUE, ellipse.level = 0.95) + theme_classic()


#plot graph with the legend on top and without legend name
pca_plot_spermiogram + theme(legend.position = "top") + theme(legend.title=element_blank())

#count
pca_plot_count<-fviz_pca_ind(res.pca, pointsize = 3,
                       pointshape = 21, fill = metadata$count,
                       repel = TRUE, # Avoid text overlapping (slow if many points)
                       label = "none", title= "Count",
                       mean.point=F,
                       addEllipses = TRUE, ellipse.level = 0.95) + theme_classic()


#count
pca_plot_conc<-fviz_pca_ind(res.pca, pointsize = 3,
                       pointshape = 21, fill = metadata$conc,
                       repel = TRUE, # Avoid text overlapping (slow if many points)
                       label = "none", title= "concentration",
                       mean.point=F,
                       addEllipses = TRUE, ellipse.level = 0.95) + theme_classic()
                       
#prog_mot
pca_plot_progmot<-fviz_pca_ind(res.pca, pointsize = 3,
                       pointshape = 21, fill = metadata$prog_mot,
                       repel = TRUE, # Avoid text overlapping (slow if many points)
                       label = "none", title= "Progressive motility",
                       mean.point=F,
                       addEllipses = TRUE, ellipse.level = 0.95) + theme_classic()
                        
#tot_mot
pca_plot_totmot<-fviz_pca_ind(res.pca, pointsize = 3,
                               pointshape = 21, fill = metadata$tot_mot,
                               repel = TRUE, 
                               label = "none", title= "Total motility",
                               mean.point=F,
                               addEllipses = TRUE, ellipse.level = 0.95) +                                          theme_classic()

#morphology
pca_plot_morph<-fviz_pca_ind(res.pca, pointsize = 3,
                              pointshape = 21, fill = metadata$morph,
                              repel = TRUE,
                              mean.point=F,
                              label = "none", title= "Total motility",
                              addEllipses = TRUE, ellipse.level = 0.95) +
                              theme_classic()
```

## Heatmap - Figure 5

```
#microbiota data
microbiota <- read.csv("~/sperm_metabolites/microbiota_data.csv")

#transform microbiota data - add row names and transpose
rownames.table<-microbiota$id
#microbiota<-subset(microbiota, select=-c(id))
microbiota$id<-NULL
row.names(microbiota)<-rownames.table
microbiota<-t(microbiota)

#Select only bacteria that are at least 1%
microbiota_above1<-colSums(microbiota)/19
microbiota_above1<-microbiota_above1[microbiota_above1<0.01]
microbiota_above1<-names(microbiota_above1)
microbiota<-microbiota[ , !(colnames(microbiota) %in% microbiota_above1)]

#add row names to metadata table
metadata <- read.csv("~/sperm_metabolites/input_metadata.csv")
rownames.table<-metadata$id
row.names(metadata)<-rownames.table
metadata$id<-NULL

#select spermiogram data from metadata
spermiogram<-subset(metadata,select=c(conc_n,count_n,prog_mot_n,tot_mot_n))

#metabolites data from the volcano plot part - add row names and transpose

#SWITCH
metabolites<-metab_count

metabolites$metab<-NULL
metabolites<-ddply(metabolites,"taxonomy",numcolwise(sum))
length.db<-dim(metabolites)[1]
if(is.na(metabolites$taxonomy[length.db])) {
 metabolites<-metabolites[-(length.db),]
 }  
row.names(metabolites)<-metabolites$taxonomy
metabolites$taxonomy<-NULL

# row.names(metabolites)<-metabolites$metab
# metabolites<-subset(metabolites, select=-c(metab))
tmetabolites<-t(metabolites)

#normalize metabolites and spermiogram data
# norm_metabolites = log10(tmetabolites+0.000001)
# norm_spermiogram = log10(spermiogram+0.000001)

norm_metabolites = tmetabolites
norm_spermiogram = log2(spermiogram)

#SWITCH  - metabolites-microbiota or metabolites-spermiogram
#correlation<-cor(microbiota,norm_metabolites, method = "kendall")
correlation<-cor(microbiota,norm_metabolites, method = "spearman")
tcorrelation<-t(correlation)


plot(sort(microbiota[,10]))

#correlation<-cor(spermiogram, microbiota, method = "spearman")
#tcorrelation<-t(correlation)

####CREATE HEATMAP WITH MOST INFORMATIVE COLUMNS 1######

new_correlation<-matrix(nrow = 16, ncol = 0)

for (n in 1:ncol(correlation)) {
  if(max(abs(correlation[,n]))>0.5){
    
    new_correlation<-cbind(new_correlation,correlation[,n])
    
    colnames(new_correlation)[ncol(new_correlation)] <- colnames(correlation)[n]
  }
}

#create heatmap
pheatmap(new_correlation, scale = "row",cellwidth = 10, cellheight = 9,fontsize = 9)
plot.background = element_rect(fill = NA, color = NA)

# Pairwise correlation between samples (columns)
cols.cor_new <- cor(new_correlation, use = "pairwise.complete.obs", method = "spearman")
# Pairwise correlation between rows (genes)
rows.cor_new <- cor(t(new_correlation), use = "pairwise.complete.obs", method = "spearman")

heatmap<-pheatmap(
  new_correlation, scale = "row", cellwidth = 8, cellheight = 8,fontsize = 8,
  clustering_distance_cols = as.dist(1 - cols.cor_new),
  clustering_distance_rows = as.dist(1 - rows.cor_new),
  background = "transparent",cutree_rows = 4)

#save plot if needed
pdf("~/sperm_metabolites/output/heatmap_microbiota_metab_spearman_all.pdf", width = 20,height = 8,bg = "transparent")
heatmap
dev.off()



#add barplot of taxa percentage following heatmap entries order
tot_microbiota<-colSums(microbiota)
tot_microbiota<-as.data.frame(tot_microbiota)
percentage_microbiota<-tot_microbiota/(colSums(tot_microbiota))*100
percentage_microbiota<-as.data.frame(percentage_microbiota)
percentage_microbiota_number<-percentage_microbiota
percentage_microbiota_number[,2]<-c(1:16)
genus_order<-heatmap$tree_row$order

#heatmap$tree_col$labels

library("tibble")
percentage_microbiota_ord <- as_tibble(percentage_microbiota)
percentage_microbiota_ord<-t(percentage_microbiota_ord)
percentage_microbiota_ord <- percentage_microbiota_ord[,genus_order]

barplot(percentage_microbiota_ord,ylim=c(30,0))


#add barplot of metabolite quantification
tot_metabolite<-rowSums(metabolites)
tot_metabolite<-data.frame(tot_metabolite)
tot_metabolite$metab<-row.names(tot_metabolite)
colnames(tot_metabolite)<-c("value","metab")

metabolite_order<-data.frame(heatmap$tree_col$labels)
colnames(metabolite_order)<-"metab"
metabolite_order<-heatmap$tree_col$labels
metabolite_order<-data.frame(metabolite_order[heatmap$tree_col$order])
colnames(metabolite_order)<-"metab"

metabolite_table <- left_join(metabolite_order,tot_metabolite,by="metab")
row.names(metabolite_table)<-metabolite_table$metab

#plot  
metab_data<-ggplot(metabolite_table,aes(x=factor(metab, level = rev(metabolite_table$metab)),value)) +
  geom_bar(stat="identity",fill="lightgrey") +
  #geom_text(aes(x = metab, y = value, label = metabolite_table$metab),hjust = -0.1,size=3)+
  scale_y_log10() +
  theme_classic() +
  coord_flip() +
  theme(axis.text.y = element_text(hjust = 0))+
  theme(rect = element_rect(fill = "transparent"))

pdf("~/sperm_metabolites/output/metab_data.pdf", width = 8,height = 12,bg = "transparent")
print(metab_data)
dev.off()

#add barplot of microbiota quantification
tot_microbiota<-colSums(microbiota)
tot_microbiota<-as.data.frame(tot_microbiota)
percentage_microbiota<-tot_microbiota/(colSums(tot_microbiota))*100
percentage_microbiota<-as.data.frame(percentage_microbiota)
percentage_microbiota_number<-percentage_microbiota
percentage_microbiota_number[,2]<-c(1:16)
percentage_microbiota_number$genus<-row.names(percentage_microbiota_number)

microbiota_order<-data.frame(heatmap$tree_row$labels)
colnames(microbiota_order)<-"genus"
microbiota_order<-heatmap$tree_row$labels
microbiota_order<-data.frame(microbiota_order[heatmap$tree_row$order])
colnames(microbiota_order)<-"genus"

microbiota_table <- left_join(microbiota_order,percentage_microbiota_number,by="genus")
row.names(microbiota_table)<-microbiota_table$genus

#plot  
microb_data<-ggplot(microbiota_table,aes(x=factor(genus, level = microbiota_table$genus),tot_microbiota)) +
  geom_bar(stat="identity",fill="lightgrey") +
  #geom_text(aes(x = metab, y = value, label = metabolite_table$metab),hjust = -0.1,size=3)+
  theme_classic() +
  coord_flip() +
  theme(axis.text.y = element_text(hjust = 0))+
  theme(rect = element_rect(fill = "transparent"))

pdf("~/sperm_metabolites/output/microb_data.pdf", width = 8,height = 10,bg = "transparent")
print(microb_data)
dev.off()

###create table with significance for the heatmap
microbiota <- read.csv("~/sperm_metabolites/microbiota_data.csv")
row.names(microbiota)<-microbiota$id
microbiota$id<-NULL

#Select only bacteria that are at least 1%
microbiota_above1<-colSums(microbiota)/19
microbiota_above1<-microbiota_above1[microbiota_above1<0.01]
microbiota_above1<-names(microbiota_above1)
microbiota<-microbiota[!(rownames(microbiota) %in% microbiota_above1),]

#SWITCH
input<-metab_count

input$metab<-NULL
input<-ddply(input,"taxonomy",numcolwise(sum))
length.db<-dim(input)[1]
input$taxonomy[length.db]<-"NA"
colnames(input)[1]<-"metab"

row.names(input)<-input$metab
input$metab<-NULL

input<-log10(input)

microbiota_corr<-t(as.matrix(microbiota))

#loop
coeff_table<-data.frame(matrix(NA, nrow = 0, ncol = 4))

for (z in 1:ncol(microbiota_corr)) {
  
  for (i in 1:nrow(input)) {
    
    test<-as.matrix(rbind(input[i,], microbiota_corr[,z]))
    test<-t(test)
    test<-as.data.frame(test)
    colnames(test)[2]<-colnames(microbiota_corr)[z]
    
    #test1<-data.frame(matrix(NA, nrow = 0, ncol = 2))
    
    test$V1<-as.numeric(as.character(unlist(test[,1])))
    test$V2<-as.numeric(as.character(unlist(test[,2])))
    
    colnames(test)[3]<-"V1"
    colnames(test)[4]<-"V2"
    
    res <- cor.test(test$V1, test$V2, 
                    method = "spearman")
    
    coeff_table<-rbind(coeff_table,c(colnames(test)[1],colnames(test)[2],res$estimate,res$p.value))
    colnames(coeff_table)<-c("metab","bacterium","coeff","p-value")
  }
}

coeff_table_p<-coeff_table[,(-3)]
colnames(coeff_table_p)<-c("metab","bacterium","pvalue")
significance<-cast(coeff_table_p, bacterium~metab)
row.names(significance)<-significance$bacterium
significance$bacterium<-NULL


significance<-significance[heatmap$tree_row$labels,heatmap$tree_col$labels]
significance_table<-replace(significance, significance < 0.05, "*")
significance_table<-replace(significance_table, significance_table > 0.05, "") 

heatmap_sign<-pheatmap(
  new_correlation, scale = "row", cellwidth = 8, cellheight = 8,fontsize = 8,
  clustering_distance_cols = as.dist(1 - cols.cor_new),
  clustering_distance_rows = as.dist(1 - rows.cor_new),
  background = "transparent",cutree_rows = 4,
  display_numbers = significance_table)

#save plot if needed
pdf("~/sperm_metabolites/output/heatmap_microbiota_metab_spearman_all_sign.pdf", width = 20,height = 8,bg = "transparent")
heatmap_sign
dev.off()

test<-test[,metabolite_table$metab]
test<-test[microbiota_table$genus,]

significance_table<-replace(test, test < 0.05, "*")
significance_table<-replace(significance_table, test > 0.05, "") 

write.csv(significance_table,"significance_table_heatmap.csv")

sign_table<-significance_table
colnames(sign_table)<-c(10:60)
rownames(sign_table)<-c(1:16)

pdf("significance_table_heatmap.pdf",width = 20,height = 12,bg = "transparent")
grid.table(sign_table)
dev.off()

heatmap$tree_row$order


#save plot if needed
pdf("~/sperm_metabolites/output/heatmap_microbiota_metab_spearman_all.pdf", width = 20,height = 8,bg = "transparent")
heatmap
dev.off()

display_numbers = test_labels
```

## Urocanate in vitro testing - Figure 6

```
viability<-read.csv("~/sperm_metabolites/urocanate_physiology.csv",header=T)

viability_table<-ddply(viability, c("sample"), summarise,
                     mean = mean(value), sd = sd(value))

viability_table<-cbind(viability_table,color=c("#619CFF","#F8766D","#619CFF","#F8766D","#619CFF","#F8766D","#619CFF","#F8766D"),group=c("A","A","B","B","C","C","D","D"))


bar_plot_viability <- ggplot(viability_table, aes(x=sample,y=mean)) +
  geom_bar(stat = "identity",fill = viability_table$color,position='dodge') +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  labs(x = "", y = "Viable spermatozoa") +
  #scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  theme_minimal() +
  coord_cartesian(ylim=c(0,1)) +
  theme_classic() +
  scale_x_discrete(labels=c("control","treated","control","treated","control","treated","control","treated")) +
  theme(text = element_text(size=12),axis.text.x = element_text(angle = 45, hjust=1)) +facet_grid(~viability_table$group,scales="free")


mmp<-read.csv("~/sperm_metabolites/urocanate_mmp.csv",header=T)

mmp_table<-ddply(mmp, c("sample"), summarise,
                     mean = mean(value), sd = sd(value))

mmp_table<-cbind(mmp_table,color=c("#619CFF","#F8766D","#619CFF","#F8766D","#619CFF","#F8766D","#619CFF","#F8766D"),group=c("A","A","B","B","C","C","D","D"))


bar_plot_mpp <- ggplot(mmp_table, aes(x=sample,y=mean)) +
  geom_bar(stat = "identity",fill = mmp_table$color,position='dodge') +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
  labs(x = "", y = "Mitochondiral membrane potential") +
  #scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  theme_minimal() +
  #coord_cartesian(ylim=c(0,1)) +
  theme_classic() +
  scale_x_discrete(labels=c("control","treated","control","treated","control","treated","control","treated")) +
  theme(text = element_text(size=12),axis.text.x = element_text(angle = 45, hjust=1))+
  facet_grid(~mmp_table$group,scales="free")

grid.arrange(bar_plot_viability,bar_plot_mpp,ncol=2)
```
