# Download data programmatically

Whenever a paper is published, genomic data almost also must be made
publicly avaialble through an online repository. There are several
commonly used ones, but perhaps the most ubiquitious one is the [Gene
Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/), or simply GEO.
Data are assigned accession number that are often found in the “data
availability” section of published papers. These accession numbers lead
you to where the processed data are stored on GEO. (As an aside, the raw
sequencing data tied to this accession number as well, but the raw data
are stored in the Sequence Read Archive, which is beyond the scope of
this tutorial.)

As an example, the accession number
[GSE139324](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139324)
takes you to the page where the processed scRNAseq data are stored for
our 2020 Immunity paper.

If you go down to the *Samples* section, and click on eg
[GSM4138110](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4138110),
it will take you to a description of that sample and links to download
the files individually. So, to look at these data in Seurat, you’d have
to:

-   Download these 3 files
    -   GSM4138110\_HNSCC\_1\_PBMC\_barcodes.tsv.gz
    -   GSM4138110\_HNSCC\_1\_PBMC\_genes.tsv.gz  
    -   GSM4138110\_HNSCC\_1\_PBMC\_matrix.mtx.gz
-   Put them into a folder
-   Read them into Seruat

Not too bad, right? It’s fine for a couple of samples, but if you want
to download all 63 samples from this study it gets tedious quickly.

Instead, we can write a function and use the Bioconductor R package
[GEOquery](https://bioconductor.org/packages/release/bioc/html/GEOquery.html).

# Prerequisites

To run this tutorial, we assume the following:

-   You have downloaded and installed the GEOquery and Biobase packages
    through Bioconductor.
-   You have some familiarity with the concept of writing functions
    -   For a refresher on writing function, [check this
        out](https://r4ds.hadley.nz/functions#introduction) from R for
        Data Science.

# Load packages

First, we will load the necessary packages.

    library(GEOquery)

    ## Loading required package: Biobase

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Setting options('download.file.method.GEOquery'='auto')

    ## Setting options('GEOquery.inmemory.gpl'=FALSE)

    library(Biobase)
    library(dplyr)

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     combine

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

# Identify samples and metadata from an accession number

Here, we will gain some familiarity with the data structures that we can
pull from GEO and how we can download samples using this information.

    gds_use <- GEOquery::getGEO("GSE139324")

    ## Found 1 file(s)

    ## GSE139324_series_matrix.txt.gz

    gds_use

    ## $GSE139324_series_matrix.txt.gz
    ## ExpressionSet (storageMode: lockedEnvironment)
    ## assayData: 0 features, 63 samples 
    ##   element names: exprs 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: GSM4138110 GSM4138111 ... GSM4138172 (63 total)
    ##   varLabels: title geo_accession ... tissue:ch1 (44 total)
    ##   varMetadata: labelDescription
    ## featureData: none
    ## experimentData: use 'experimentData(object)'
    ##   pubMedIds: 31924475
    ## 34099645 
    ## Annotation: GPL18573

    class(gds_use) # list

    ## [1] "list"

    class(gds_use[[1]]) # first item in the list 

    ## [1] "ExpressionSet"
    ## attr(,"package")
    ## [1] "Biobase"

    slotNames(gds_use[[1]]) # all the containers for the data in this S4 object

    ## [1] "experimentData"    "assayData"         "phenoData"        
    ## [4] "featureData"       "annotation"        "protocolData"     
    ## [7] ".__classVersion__"

Let’s unpack what we’re seeing in the ExpressionSet object.

We have 7 slots that ostensibly are holding some sorts of data. But we
can see that assayData has 0 features so is probably not helpful. We can
also see that protocolData has nothing in it and neither does
featureData.

There are some other things here, but it looks like all the goodies are
in the phenoData slot.

# Explore the phenoData

    gds_use[[1]]@phenoData # using the "@" to access slots in the ExpressionSet

    ## An object of class 'AnnotatedDataFrame'
    ##   sampleNames: GSM4138110 GSM4138111 ... GSM4138172 (63 total)
    ##   varLabels: title geo_accession ... tissue:ch1 (44 total)
    ##   varMetadata: labelDescription

    class(gds_use[[1]]@phenoData) # AnnotatedDataFrame

    ## [1] "AnnotatedDataFrame"
    ## attr(,"package")
    ## [1] "Biobase"

    # Looks like individual sample accessions are under sampleNames
    # varLabels has 44 different pieces of data

    gds_use[[1]]@phenoData@data %>%
      glimpse()

    ## Rows: 63
    ## Columns: 44
    ## $ title                     <chr> "HNSCC_1_PBMC", "HNSCC_1_TIL", "HNSCC_2_PBMC…
    ## $ geo_accession             <chr> "GSM4138110", "GSM4138111", "GSM4138112", "G…
    ## $ status                    <chr> "Public on Nov 13 2019", "Public on Nov 13 2…
    ## $ submission_date           <chr> "Oct 24 2019", "Oct 24 2019", "Oct 24 2019",…
    ## $ last_update_date          <chr> "Nov 13 2019", "Nov 13 2019", "Nov 13 2019",…
    ## $ type                      <chr> "SRA", "SRA", "SRA", "SRA", "SRA", "SRA", "S…
    ## $ channel_count             <chr> "1", "1", "1", "1", "1", "1", "1", "1", "1",…
    ## $ source_name_ch1           <chr> "peripheral blood", "tumor tissue", "periphe…
    ## $ organism_ch1              <chr> "Homo sapiens", "Homo sapiens", "Homo sapien…
    ## $ characteristics_ch1       <chr> "tissue: peripheral blood", "tissue: tumor t…
    ## $ characteristics_ch1.1     <chr> "disease state: Head and neck squamous cell …
    ## $ characteristics_ch1.2     <chr> "hpv_status: HPV_negative", "hpv_status: HPV…
    ## $ treatment_protocol_ch1    <chr> "Ficoll density gradient centrifugation to i…
    ## $ molecule_ch1              <chr> "polyA RNA", "polyA RNA", "polyA RNA", "poly…
    ## $ extract_protocol_ch1      <chr> "10x 3' v2 scRNAseq", "10x 3' v2 scRNAseq", …
    ## $ taxid_ch1                 <chr> "9606", "9606", "9606", "9606", "9606", "960…
    ## $ data_processing           <chr> "Demultiplexing using CellRanger mkfastq", "…
    ## $ data_processing.1         <chr> "Alignment and generation of gene/barcode ma…
    ## $ data_processing.2         <chr> "Genome_build: GRCh38", "Genome_build: GRCh3…
    ## $ data_processing.3         <chr> "Supplementary_files_format_and_content: fil…
    ## $ platform_id               <chr> "GPL18573", "GPL18573", "GPL18573", "GPL1857…
    ## $ contact_name              <chr> "Anthony,Richard,Cillo", "Anthony,Richard,Ci…
    ## $ contact_email             <chr> "arc85@pitt.edu", "arc85@pitt.edu", "arc85@p…
    ## $ contact_laboratory        <chr> "Vignali Lab", "Vignali Lab", "Vignali Lab",…
    ## $ contact_department        <chr> "Immunology", "Immunology", "Immunology", "I…
    ## $ contact_institute         <chr> "University of Pittsburgh", "University of P…
    ## $ contact_address           <chr> "Suite 2.19 Hillman Cancer, 5115 Centre Ave"…
    ## $ contact_city              <chr> "Pittsburgh", "Pittsburgh", "Pittsburgh", "P…
    ## $ contact_state             <chr> "PA", "PA", "PA", "PA", "PA", "PA", "PA", "P…
    ## $ `contact_zip/postal_code` <chr> "15232", "15232", "15232", "15232", "15232",…
    ## $ contact_country           <chr> "USA", "USA", "USA", "USA", "USA", "USA", "U…
    ## $ data_row_count            <chr> "0", "0", "0", "0", "0", "0", "0", "0", "0",…
    ## $ instrument_model          <chr> "Illumina NextSeq 500", "Illumina NextSeq 50…
    ## $ library_selection         <chr> "cDNA", "cDNA", "cDNA", "cDNA", "cDNA", "cDN…
    ## $ library_source            <chr> "transcriptomic", "transcriptomic", "transcr…
    ## $ library_strategy          <chr> "RNA-Seq", "RNA-Seq", "RNA-Seq", "RNA-Seq", …
    ## $ relation                  <chr> "BioSample: https://www.ncbi.nlm.nih.gov/bio…
    ## $ relation.1                <chr> "SRA: https://www.ncbi.nlm.nih.gov/sra?term=…
    ## $ supplementary_file_1      <chr> "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4…
    ## $ supplementary_file_2      <chr> "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4…
    ## $ supplementary_file_3      <chr> "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4…
    ## $ `disease state:ch1`       <chr> "Head and neck squamous cell carcinoma (HNSC…
    ## $ `hpv_status:ch1`          <chr> "HPV_negative", "HPV_negative", "HPV_negativ…
    ## $ `tissue:ch1`              <chr> "peripheral blood", "tumor tissue", "periphe…

There is a lot of info here, but what eventually becomes clear is the
following: - There’s lots of useful metadata here for each sample -
There are ftp links to supplementary files (in this case, the 3 we need
for scRNAseq) that we can download

# Write a function to download data from GEO

We will now build a function that uses the info we can get from the
ExpressionSet phenoData to download a few samples from this study.

    geo_download <- function(target_directory,geo_accession,samples_to_download) {

        # Call query to GEO
        gds <- GEOquery::getGEO(geo_accession)

        # Check whether the files are already present in the target_directory
        file_names <- list.files(target_directory)
        num_samples_present <- length(which(samples_to_download %in% file_names))

        if (num_samples_present==length(samples_to_download)) {

            print(
              "All samples already present!"
            )

        } else if (num_samples_present>0) {
          
           print(paste(
                "There are ",
                num_samples_present,
                " of ",
                length(samples_to_download),
                " already here!",
                sep=""
            )
            )
          
          sub_samples_to_download <- setdiff(samples_to_download,file_names)

          # Identify samples from GEO
          data_use <- Biobase::phenoData(gds[[1]])@data[Biobase::phenoData(gds[[1]])@data$title %in% sub_samples_to_download,]

            # Create subdirectories and download from GEO
            for (i in 1:nrow(data_use)) {
        
                dir_name <- paste(target_directory,data_use$title[i],sep="/")
                dir.create(dir_name)

                tmp_dwnload1 <- as.character(data_use[i,which(grepl("supplementary_file",colnames(data_use)))[1]])
                tmp_dwnload2 <- as.character(data_use[i,which(grepl("supplementary_file",colnames(data_use)))[2]])
                tmp_dwnload3 <- as.character(data_use[i,which(grepl("supplementary_file",colnames(data_use)))[3]])

                download.file(tmp_dwnload1,cacheOK = FALSE,destfile=paste(dir_name,"barcodes.tsv.gz",sep="/"))
                download.file(tmp_dwnload2,cacheOK = FALSE,destfile=paste(dir_name,"features.tsv.gz",sep="/"))
                download.file(tmp_dwnload3,cacheOK = FALSE,destfile=paste(dir_name,"matrix.mtx.gz",sep="/"))
            }

        } else {

            # Identify samples from GEO
            data_use <- Biobase::phenoData(gds[[1]])@data[Biobase::phenoData(gds[[1]])@data$title %in% samples_to_download,]

            # Create subdirectories and download from GEO
            for (i in 1:nrow(data_use)) {
        
                dir_name <- paste(target_directory,data_use$title[i],sep="/")
                dir.create(dir_name)

                tmp_dwnload1 <- as.character(data_use[i,which(grepl("supplementary_file",colnames(data_use)))[1]])
                tmp_dwnload2 <- as.character(data_use[i,which(grepl("supplementary_file",colnames(data_use)))[2]])
                tmp_dwnload3 <- as.character(data_use[i,which(grepl("supplementary_file",colnames(data_use)))[3]])

                download.file(tmp_dwnload1,cacheOK = FALSE,destfile=paste(dir_name,"barcodes.tsv.gz",sep="/"))
                download.file(tmp_dwnload2,cacheOK = FALSE,destfile=paste(dir_name,"features.tsv.gz",sep="/"))
                download.file(tmp_dwnload3,cacheOK = FALSE,destfile=paste(dir_name,"matrix.mtx.gz",sep="/"))
            }
          
        }

    }

## Try out the function

Running this function for the first time will download the data in the
proper structure for it to be read directly into Seurat.

Attempting to run the function a section time will tell you that you’ve
already downloaded the data!

    geo_download(target_directory="~/Desktop/geo_download_test",
        geo_accession="GSE139324",
        samples_to_download=c("HD_PBMC_1","HD_PBMC_2","HD_PBMC_3")
        )

    ## Found 1 file(s)

    ## GSE139324_series_matrix.txt.gz

    ## Using locally cached version: /var/folders/1w/32zp8ffj6q50nny8nt6y1bvc0000gq/T//RtmpBGdiNR/GSE139324_series_matrix.txt.gz

    ## Using locally cached version of GPL18573 found here:
    ## /var/folders/1w/32zp8ffj6q50nny8nt6y1bvc0000gq/T//RtmpBGdiNR/GPL18573.soft.gz

    ## [1] "All samples already present!"

    # Will tell you all samples are already here 
    geo_download(target_directory="~/Desktop/geo_download_test",
        geo_accession="GSE139324",
        samples_to_download=c("HD_PBMC_1","HD_PBMC_2","HD_PBMC_3")
        )

    ## Found 1 file(s)

    ## GSE139324_series_matrix.txt.gz

    ## Using locally cached version: /var/folders/1w/32zp8ffj6q50nny8nt6y1bvc0000gq/T//RtmpBGdiNR/GSE139324_series_matrix.txt.gz

    ## Using locally cached version of GPL18573 found here:
    ## /var/folders/1w/32zp8ffj6q50nny8nt6y1bvc0000gq/T//RtmpBGdiNR/GPL18573.soft.gz

    ## [1] "All samples already present!"

    # Will tell you 3 of 4 are already here
    # And will download the missing file 
    geo_download(target_directory="~/Desktop/geo_download_test",
        geo_accession="GSE139324",
        samples_to_download=c("HD_PBMC_1","HD_PBMC_2","HD_PBMC_3","HD_PBMC_4")
        )

    ## Found 1 file(s)

    ## GSE139324_series_matrix.txt.gz

    ## Using locally cached version: /var/folders/1w/32zp8ffj6q50nny8nt6y1bvc0000gq/T//RtmpBGdiNR/GSE139324_series_matrix.txt.gz

    ## Using locally cached version of GPL18573 found here:
    ## /var/folders/1w/32zp8ffj6q50nny8nt6y1bvc0000gq/T//RtmpBGdiNR/GPL18573.soft.gz

    ## [1] "All samples already present!"

## Save the metadata

Let’s also save the metadata. It will be helpful for including some of
these pieces in Seurat when/if we analyze these samples.

    metadata_to_save <- gds_use[[1]]@phenoData@data %>%
      as_tibble(.,rownames="sample_accession")

    readr::write_tsv(metadata_to_save,
                     file="~/Desktop/geo_download_test/GSE139324_metadata.tsv")

## Conclusions

Here, we’ve gained some familiarity with programmatically downloading
data from GEO. This will be helpful for analysis of publicly available
data!
