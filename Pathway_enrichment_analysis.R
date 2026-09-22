suppressPackageStartupMessages({

    library(Seurat)
    library(Matrix)
    library(svglite)
    library(tidyverse)

})

##############################################################################
## PARAMETERS
##############################################################################

cluster_name <- "Cluster_G"

min_pct <- 0.005

obj_file <- "results/objects/PGintegrated73.reclustered.linkpeaks.trajectory.rds"

module_file <- file.path(
    "results",
    "temporal_modules",
    cluster_name,
    "RNA_gene_modules.csv"
)

out_dir <- file.path(
    "results",
    "temporal_ontology",
    cluster_name
)

dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

##############################################################################
## CHECK INPUTS
##############################################################################

stopifnot(file.exists(obj_file))
stopifnot(file.exists(module_file))

##############################################################################
## LOAD SEURAT OBJECT
##############################################################################

cat("\nLoading Seurat object...\n")

integrated73 <- readRDS(obj_file)

DefaultAssay(integrated73) <- "SoupXRNA"

##############################################################################
## SUBSET CLUSTER
##############################################################################

cat("Subsetting", cluster_name, "...\n")

obj <- subset(
    integrated73,
    subset = traj_cluster == cluster_name
)

cat(
    "Cells:",
    ncol(obj),
    "\n"
)

stopifnot(ncol(obj) > 0)

##############################################################################
## EXPRESSION MATRIX
##############################################################################

counts <- GetAssayData(
    obj,
    layer = "counts"
)

##############################################################################
## GENE UNIVERSE
##############################################################################

pct_expr <- Matrix::rowSums(counts > 0) / ncol(counts)

universe_genes <- rownames(counts)[
    pct_expr >= min_pct
]

cat(
    "Universe genes:",
    length(universe_genes),
    "\n"
)

stopifnot(length(universe_genes) > 100)

##############################################################################
## READ MODULE TABLE
##############################################################################

modules <- read_csv(
    module_file,
    show_col_types = FALSE
)

stopifnot(
    all(c("gene", "module") %in% colnames(modules))
)

##############################################################################
## CLEAN
##############################################################################

modules <- modules %>%
    mutate(
        gene = stringr::str_trim(gene),
        module = as.integer(module)
    ) %>%
    filter(
        !is.na(gene),
        gene != ""
    ) %>%
    distinct()

##############################################################################
## GENES BY MODULE
##############################################################################

genes_by_module <- split(
    modules$gene,
    modules$module
)

cat("\nModules detected:\n")

for(i in names(genes_by_module)){

    cat(
        "Module",
        i,
        ":",
        length(genes_by_module[[i]]),
        "genes\n"
    )

}

##############################################################################
## KEEP ONLY GENES PRESENT IN UNIVERSE
##############################################################################

genes_by_module <- lapply(
    genes_by_module,
    intersect,
    universe_genes
)

cat("\nGenes retained after universe filtering:\n")

for(i in names(genes_by_module)){

    cat(
        "Module",
        i,
        ":",
        length(genes_by_module[[i]]),
        "\n"
    )

}

stopifnot(
    all(lengths(genes_by_module) >= 5)
)




##############################################################################
## SUMMARY
##############################################################################

cat("\n----------------------------------------\n")

cat(
    "Cluster:",
    cluster_name,
    "\n"
)

cat(
    "Cells:",
    ncol(obj),
    "\n"
)

cat(
    "Universe:",
    length(universe_genes),
    "genes\n"
)

cat(
    "Modules:",
    length(genes_by_module),
    "\n"
)

cat("----------------------------------------\n")


##############################################################################
## LIBRARIES
##############################################################################

library(clusterProfiler)
library(org.Rn.eg.db)
library(org.Hs.eg.db)
library(gprofiler2)

##############################################################################
## RAT SYMBOL -> RAT ENTREZ
##############################################################################

rat_symbol_to_entrez <- function(genes){

    genes <- unique(genes)
    genes <- genes[!is.na(genes)]
    genes <- genes[genes != ""]

    if(length(genes) == 0)
        return(character(0))

    conv <- suppressMessages(

        bitr(
            genes,
            fromType = "SYMBOL",
            toType   = "ENTREZID",
            OrgDb    = org.Rn.eg.db
        )

    )

    if(is.null(conv))
        return(character(0))

    unique(conv$ENTREZID)

}

##############################################################################
## RAT SYMBOL -> HUMAN ENTREZ
##############################################################################

rat_symbol_to_human_entrez <- function(genes){

    genes <- unique(genes)
    genes <- genes[!is.na(genes)]
    genes <- genes[genes != ""]

    if(length(genes) == 0)
        return(character(0))

    orth <- tryCatch(

        gorth(

            query = genes,

            source_organism = "rnorvegicus",

            target_organism = "hsapiens"

        ),

        error = function(e) NULL

    )

    if(is.null(orth))
        return(character(0))

    if(nrow(orth) == 0)
        return(character(0))

    human_symbols <- unique(orth$ortholog_name)

    conv <- suppressMessages(

        bitr(

            human_symbols,

            fromType = "SYMBOL",

            toType = "ENTREZID",

            OrgDb = org.Hs.eg.db

        )

    )

    if(is.null(conv))
        return(character(0))

    unique(conv$ENTREZID)

}

##############################################################################
## CONVERT EACH MODULE
##############################################################################

genes_rat <- lapply(

    genes_by_module,

    rat_symbol_to_entrez

)

genes_human <- lapply(

    genes_by_module,

    rat_symbol_to_human_entrez

)

##############################################################################
## CONVERT UNIVERSE
##############################################################################

universe_rat <- rat_symbol_to_entrez(

    universe_genes

)

universe_human <- rat_symbol_to_human_entrez(

    universe_genes

)

##############################################################################
## SUMMARY
##############################################################################

cat("\n========== RAT ==========\n")

print(

    sapply(

        genes_rat,

        length

    )

)

cat(

    "\nUniverse rat:",

    length(universe_rat),

    "\n"

)

cat("\n========== HUMAN ==========\n")

print(

    sapply(

        genes_human,

        length

    )

)

cat(

    "\nUniverse human:",

    length(universe_human),

    "\n"

)

##############################################################################
## SAFETY
##############################################################################

stopifnot(

    length(universe_rat) > 5000

)

stopifnot(

    length(universe_human) > 5000

)

stopifnot(

    all(

        lengths(genes_rat) > 20

    )

)

stopifnot(

    all(

        lengths(genes_human) > 20

    )

)



##############################################################################
## ENRICHISSEMENT PER MODULE
##############################################################################

library(ReactomePA)
library(msigdbr)

##############################################################################
## HALLMARK DATABASE
##############################################################################

hallmark <- msigdbr(
    species = "Homo sapiens",
    collection = "H"
) |>
    dplyr::select(gs_name, ncbi_gene)

##############################################################################
## FONCTION
##############################################################################

run_module_enrichment <- function(
    module_name,
    rat_genes,
    human_genes
){

    cat("\n========================================\n")
    cat("Module", module_name, "\n")
    cat("Rat genes   :", length(rat_genes), "\n")
    cat("Human genes :", length(human_genes), "\n")

    ego <- tryCatch(

        enrichGO(
    gene          = human_genes,
    universe      = universe_human,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    minGSSize     = 5,
    maxGSSize     = 5000
),

        error=function(e) NULL

    )

    ekegg <- tryCatch(

 enrichKEGG(
    gene          = human_genes,
    universe      = universe_human,
    organism      = "hsa",
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    minGSSize     = 5,
    maxGSSize     = 5000
),

        error=function(e) NULL

    )

    ereact <- tryCatch(

 enrichPathway(
    gene          = human_genes,
    universe      = universe_human,
    organism      = "human",
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    minGSSize     = 5,
    maxGSSize     = 5000
),

        error=function(e) NULL

    )

    ehall <- tryCatch(

 enricher(
    gene          = human_genes,
    universe      = universe_human,
    TERM2GENE     = hallmark,
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    minGSSize     = 5,
    maxGSSize     = 5000
),

        error=function(e) NULL

    )

    list(

        GO = ego,

        KEGG = ekegg,

        Reactome = ereact,

        Hallmark = ehall

    )

}



##############################################################################
## SIMPLIFY GO
##############################################################################

simplify_all_GO <- function(
    ontology_results,
    n_keep = 200
){

    for(module in names(ontology_results)){

        go <- ontology_results[[module]]$GO

        if(is.null(go))
            next

        go_df <- as.data.frame(go)

        if(nrow(go_df) == 0)
            next

        before <- nrow(go_df)


        go_df <- go_df |>

            dplyr::arrange(p.adjust, pvalue) |>

            dplyr::slice_head(n = min(n_keep, nrow(go_df)))


        go@result <- go_df


        cat(
            "\nSimplifying GO - Module",
            module,
            "(",
            before,
            "->",
            nrow(go_df),
            "terms )..."
        )

        go <- simplify(

            go,

            cutoff = 0.5,

            by = "p.adjust",

            select_fun = min,

            measure = "Wang"

        )

        after <- nrow(as.data.frame(go))

        cat(
            " final:",
            after,
            "\n"
        )

        ontology_results[[module]]$GO <- go

    }

    ontology_results

}


ontology_results <- list()

for(module_name in names(genes_human)){

    ontology_results[[module_name]] <-

        run_module_enrichment(

            module_name,

            genes_rat[[module_name]],

            genes_human[[module_name]]

        )

}


ontology_results <- simplify_all_GO(

    ontology_results,

    n_keep = 500

)

cat("\nGO terms after simplification:\n")

print(

    sapply(

        ontology_results,

        function(x)
            nrow(as.data.frame(x$GO))

    )

)

convert_result <- function(obj, source, module){

    if(is.null(obj))
        return(NULL)

    df <- as.data.frame(obj)

    if(nrow(df)==0)
        return(NULL)

    df$Source <- source
    df$Module <- module

    df

}


enrich_tbl <- bind_rows(

lapply(

    names(ontology_results),

    function(module){

        bind_rows(

            convert_result(
                ontology_results[[module]]$GO,
                "GO",
                module
            ),

            convert_result(
                ontology_results[[module]]$KEGG,
                "KEGG",
                module
            ),

            convert_result(
                ontology_results[[module]]$Reactome,
                "Reactome",
                module
            ),

            convert_result(
                ontology_results[[module]]$Hallmark,
                "Hallmark",
                module
            )

        )

    }

)

)



enrich_tbl <- enrich_tbl |>

    dplyr::mutate(

        score = -log10(p.adjust)

    )

enrich_tbl <- enrich_tbl |>

    dplyr::mutate(

        Description = dplyr::case_when(

            Source == "Hallmark" &
                is.na(Description) ~ ID,

            TRUE ~ Description

        )

    ) |>

    dplyr::mutate(

        Description = dplyr::case_when(

            Source == "Hallmark" ~

                Description |>

                stringr::str_remove("^HALLMARK_") |>

                stringr::str_replace_all("_"," ") |>

                stringr::str_to_title(),

            TRUE ~

                Description

        )

    )


primary_tbl <- enrich_tbl |>

    dplyr::group_by(

        Source,
        ID

    ) |>

    dplyr::slice_min(

        order_by = p.adjust,

        n = 1,

        with_ties = FALSE

    ) |>

    dplyr::ungroup()

plot_terms <- primary_tbl |>

    dplyr::group_by(Module) |>

    dplyr::slice_max(

        order_by = score,

        n = 15,

        with_ties = FALSE

    ) |>

    dplyr::ungroup() |>

    dplyr::rename(

        PrimaryModule = Module,

        PrimaryScore = score,

        PrimaryQ = p.adjust

    )

plot_tbl <- enrich_tbl |>

    dplyr::inner_join(

        plot_terms |>

            dplyr::select(

                Source,
                ID,
                PrimaryModule,
                PrimaryScore,
                PrimaryQ

            ),

        by = c(

            "Source",
            "ID"

        )

    )


plot_tbl <- plot_tbl |>

    dplyr::mutate(

        GeneRatioNum = purrr::map_dbl(

            GeneRatio,

            function(x){

                y <- strsplit(x,"/")[[1]]

                as.numeric(y[1]) /
                as.numeric(y[2])

            }

        )

    )


term_order <- plot_terms |>

    dplyr::arrange(

        as.numeric(PrimaryModule),

        dplyr::desc(PrimaryScore)

    ) |>

    dplyr::pull(

        Description

    )

plot_tbl$Description <- factor(

    plot_tbl$Description,

    levels = rev(unique(term_order))

)


module_levels <- sort(
    unique(
        as.character(plot_tbl$Module)
    )
)

plot_tbl$Module <- factor(
    plot_tbl$Module,
    levels = module_levels
)

plot_tbl$PrimaryModule <- factor(
    plot_tbl$PrimaryModule,
    levels = module_levels
)

plot_tbl$Source <- factor(

    plot_tbl$Source,

    levels = c(

        "Hallmark",
        "KEGG",
        "Reactome",
        "GO"

    )

)


cat("\n========================================\n")

cat("Total enrichments :",nrow(enrich_tbl),"\n")

cat("Primary terms :",nrow(plot_terms),"\n")

cat("Dotplot points :",nrow(plot_tbl),"\n")

cat("\nPrimary terms per module\n")

print(

    table(

        plot_terms$PrimaryModule

    )

)

cat("\nDotplot points per module\n")

print(

    table(

        plot_tbl$Module

    )

)

cat("\nDotplot points by source\n")

print(

    table(

        plot_tbl$Source

    )

)


dup <- plot_terms |>

    dplyr::count(

        Description

    ) |>

    dplyr::filter(

        n > 1

    )

cat("\nDuplicated descriptions : ",nrow(dup),"\n")


write.csv(

    enrich_tbl,

    file.path(

        out_dir,

        "ontology_results_full.csv"

    ),

    row.names = FALSE

)

write.csv(

    plot_tbl,

    file.path(

        out_dir,

        "ontology_results.csv"

    ),

    row.names = FALSE

)




##############################################################################
## DOTPLOT
##############################################################################

library(ggplot2)

p <- ggplot(

    plot_tbl,

    aes(

        x = Module,

        y = Description,

        colour = score,

        size = GeneRatioNum

    )

) +

geom_point(

    alpha = 1

) +

scale_x_discrete(
    labels = paste("Module", module_levels)
) +

scale_colour_gradientn(

    colours = c(

        "#3B4CC0",
        "#7B3294",
        "#D73027"

    ),

    values = scales::rescale(

        c(

            min(plot_tbl$score),

            median(plot_tbl$score),

            max(plot_tbl$score)

        )

    ),

    name = expression(-log[10](FDR))

) +

scale_size_continuous(

    range = c(1.5,7),

    name = "Gene ratio"

) +

labs(

    x = NULL,

    y = NULL

) +

theme_bw(

    base_size = 10

) +

theme(

    panel.grid.major.x = element_blank(),

panel.grid.major.y = element_line(

    colour = "grey90",

    linewidth = 0.35

),

    axis.text.y = element_text(

        size = 7

    ),

    axis.text.x = element_text(

        size = 11,

        face = "bold"

    ),

    axis.ticks.y = element_blank(),

    axis.ticks.x = element_blank(),

    legend.title = element_text(

        size = 10,

        face = "bold"

    ),

    legend.text = element_text(

        size = 9

    ),

    plot.margin = margin(

        5,

        15,

        5,

        25

    )

)

print(p)

##############################################################################
## EXPORT
##############################################################################

ggsave(

    file.path(

        out_dir,

        "Cluster_G_temporal_ontology_dotplot.svg"

    ),

    p,

    width = 11,

    height = 15,

    limitsize = FALSE

)

ggsave(

    file.path(

        out_dir,

        "Cluster_G_temporal_ontology_dotplot.png"

    ),

    p,

    width = 11,

    height = 15,

    dpi = 600,

    limitsize = FALSE

)

ggsave(

    file.path(

        out_dir,

        "Cluster_G_temporal_ontology_dotplot.tiff"

    ),

    p,

    width = 11,

    height = 15,

    dpi = 600,

    compression = "lzw",

    limitsize = FALSE

)