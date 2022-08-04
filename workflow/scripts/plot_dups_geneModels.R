args = commandArgs(trailingOnly=TRUE)
flag_labels = list(RUN_NAME = '--run_name' , 
                   GENE_BED = '--gene_bed' ,
                   DUP_BED = "--duplicon_bed",
                   help = "-h",
                   GENE_PATTERN = "--gene_pattern",
                   OUT_FILE = "--outfile_name",
                   RLIB_PATH = "--rlib" ,
                   GENE_ANNOTATIONS = "--gene_bed_annotations"
)

################
##Process help request
################
help_stats = grep( as.character(flag_labels[['help']]) , args)
if( length(help_stats == 1) ){
  print("flags to pass:")
  print("--run_name : REQUIRED :  Run name of a given run from minigraph.smk.")
  print("--outfile_name : REQUIRED : name of pdf file you wish to output. Please include .pdf extension in name.")
  print("--gene_bed : Optional :  bed of gene models to include as a track in squashed dot plot.")
  print("--gene_bed_annotations : Optional : if you want to color genes in a specific way: pass table with gene name and some factoring variable in second column")
  print("--duplicon_bed : Optional : bed of duplicons created by dupmasker to include as track in squashed dot plot.")
  print('--gene_pattern : Optional: if --gene_bed contains all gene annotations and you want to focus on a specific gene, specify that pattern to search for here.')
  print("--rlib : Optional : if there's a R library you want to include in the path to access required packages, add it here.")
  quit()
}

#'get variables passed with flags to command line.
get_flag_var = function(args, flag_name){
  i = grep(flag_name , args)
  if( length(i) > 0 ){
    if( "=" %in% args[i] ){
      return(args[i])
    }else{
      tryCatch(
        { 
          return(args[i + 1]) 
        },
        error=function(e){
          message("Wasn't able to argument linked to flag: {flag_name} Expected whitespace between. Maybe that's the problem.")
          message(e)
          return(NA)
        },finally = {
          #do nothing
        }
      )
    }
  }else{
    return(NA)
  }
}
################
##Get flag info
################
RUN_NAME = get_flag_var(args, flag_name = flag_labels[['RUN_NAME']] )
GENE_BED = get_flag_var(args, flag_name = flag_labels[['GENE_BED']] )
DUP_BED = get_flag_var(args, flag_name = flag_labels[['DUP_BED']] )
GENE_ANNOTATIONS = get_flag_var(args, flag_name = flag_labels[['GENE_ANNOTATIONS']] )
GENE_PATTERN = get_flag_var(args, flag_name = flag_labels[['GENE_PATTERN']] )
OUT_FILE = get_flag_var(args, flag_name = flag_labels[['OUT_FILE']] )
RLIB_PATH = get_flag_var(args, flag_name = flag_labels[['RLIB_PATH']] )


################
##load libararies
################
if(!is.na(RLIB_PATH)){
  .libPaths( c( .libPaths(), RLIB_PATH ) )
}
if(! require("R.utils")) install.packages("R.utils")
if(! require("tidyverse")) install.packages("tidyverse")
if(! require("data.table")) install.packages("data.table")
if(! require("glue")) install.packages("glue")
if(! require("RColorBrewer")) install.packages("RColorBrewer")
if(! require("scales")) install.packages("scales")
if(! require("cowplot")) install.packages("cowplot")
if(! require("gridExtra")) install.packages("gridExtra")



#'tri_bed : generate data.table for geom_polygon triangles from bed_file
#'scrapped from Mitchell's Minigraph.R pipeline, this will allow me to generate polygons for duplicons and gene models.
#'@param: f : fstring of bed file to be processed
#'@param: s : size of triangles : default: 0.2
#'@param: y_mid : midpoint of polygon: important for having various tracks on a plot and making sure they don't overlap
tri_bed <- function(f, s=.2, y_mid = 0.5){
  #df = fread(glue("../sd_regions_in_hifi_wga/lpa/minimiro/temp_minimiro/{gene}_query.fasta.duplicons.extra")); df
  df=fread(f)
  names(df)[1:3]=c("query","start","end")
  if(length(names(df)) >= 6){
    names(df)[6] = "strand"
  }else{
    df$strand = "+"
  }
  df = df[order(query,start)]
  df$query = factor(df$query)
  df$y = y_mid
  df$tri_id = 1:nrow(df)
  zs = s/2 #size of smaller end of triangle
  data.table(df %>% 
               rowwise() %>% 
               mutate(xs=list(c(start, start, end, end)), 
                      ys = case_when(
                        strand == "+" ~ list(c(y+s,y-s,y-zs,y+zs)),
                        strand == "-" ~ list(c(y-zs,y+zs,y+s,y-s))
                      )) %>% unnest(cols=c("xs","ys")))
}


DUPLICONS = tri_bed( DUP_BED , s = 0.3, y_mid = 1.3)
DUPLICONS = DUPLICONS %>% 
  mutate( run = RUN_NAME, query = sub(pattern = glue("{RUN_NAME}_*"), replacement = "", x = query) ) %>% 
  separate( col = query, sep = "__", into = c("q_samp", "q_hap", NA)) %>%
  relocate(q_samp, .before = start) %>% relocate(q_hap, .after = q_samp ) %>% 
  mutate(query = paste(q_samp, q_hap, sep = "__"))
#change RGB to hex color
DUPLICONS$color = sapply(DUPLICONS$color, FUN = function(cur_rgb){ 
  cur_rgb = as.numeric( unlist(strsplit(cur_rgb, split = ",")) )
  rgb(cur_rgb[1], cur_rgb[2], cur_rgb[3], maxColorValue = 255)
} )


genes = tibble(read_tsv(GENE_BED), s = 0.3, y_mid = 0.5)
ALL_GENES = tri_bed( GENE_BED , s = 0.3, y_mid = 0.5)
names(ALL_GENES)[4] = "gene_name"
#add run_name ; process query name ; add size ; 
ALL_GENES = ALL_GENES %>% 
  mutate( run = RUN_NAME , query = sub(pattern = glue("{RUN_NAME}_*"), replacement = "" , x = query) ) %>% 
  separate( col = query , sep = "__" , into = c("q_samp" , "q_hap" ,NA )) %>% 
  mutate(query = paste(.$q_samp, .$q_hap , sep = "__")) %>%
  mutate(size = end - start) %>%
  relocate( query , .before = q_samp )
# filter out different isoforms of same gene ; take top 4 (b/c polygon)
if(!is.na(GENE_PATTERN)){ ##Process A GENE_PATTERN --> gff with many genes
  ALL_GENES = ALL_GENES %>% filter( grepl( gene_name ,pattern = GENE_PATTERN) )
  ALL_GENES = ALL_GENES %>% group_by(q_samp, q_hap , gene_name) %>% filter(size == max(size)) %>% slice_head(n = 4) %>% ungroup()
}

#add coloring to genes
gene_annotations = tibble(read_tsv(GENE_ANNOTATIONS))
ALL_GENES = left_join(ALL_GENES, gene_annotations, by = c("gene_name" = "dup"))
#figure out colors:
length(unique(ALL_GENES$cluster))


plt = ggplot(data = DUPLICONS) + 
  geom_polygon( aes(x=xs, y=ys, group=tri_id), 
                fill = DUPLICONS$color) + 
  facet_wrap(~query , ncol = 1)

plt = plt + geom_polygon(data=ALL_GENES, 
                         aes(x=xs, y=ys, group=tri_id, fill = as.factor(cluster) ) ) + # , fill = "red") + 
  facet_wrap(~query , ncol = 1) + 
  scale_fill_brewer(palette="Spectral")
  

n_samples = length( unique(DUPLICONS$query) ) 
pdf_height = n_samples * 1
pdf_width = 10
pdf(OUT_FILE, onefile = TRUE, height = pdf_height, width = pdf_width)
plt
dev.off()

