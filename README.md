# scRNAseq

## Gene-level summaries of public single-cell RNA-seq datasets

This package contains a collection of four publicly available single-cell RNA-seq data. The data were downloaded from NCBI's SRA or from EBI's ArrayExpress (see below for Accession numbers) and pre-processed with a standardized pipeline (see vignette). For each dataset, the package provides gene-level read counts obtained from Tophat + featureCounts and gene-level FPKM obtained from Cufflinks.

The dataset _fluidigm_ contains 65 cells from Pollen et al. (2014), each sequenced at high and low coverage (SRA: SRP041736).

The dataset _th2_ contains 96 T helper cells from Mahata et al. (2014) (ArrayExpress: E-MTAB-2512).

The dataset _shalek_dc_ contains 2397 dendritic cells from Shalek et al. (2014) (SRA: SRP027537).

The dataset _allen_ contains 379 cells from the mouse visual cortex. This is a subset of the data published in Tasic et al. (2016) (SRA: SRP061902).

## To install

You can download the latest release [here](https://github.com/drisso/scRNAseq/releases). Once dowloaded, type in R
```{r}
install.packages("scRNAseq_0.0.8.tar.gz", repos=NULL)
```
to install the package, where `scRNAseq_0.0.8.tar.gz` is the name of the file downloaded (note that the name could change reflecting the release of a new version of the package).

## References

Pollen, Nowakowski, Shuga, Wang, Leyrat, Lui, Li, Szpankowski, Fowler, Chen, Ramalingam, Sun, Thu, Norris, Lebofsky, Toppani, Kemp II, Wong, Clerkson, Jones, Wu, Knutsson, Alvarado, Wang, Weaver, May, Jones, Unger, Kriegstein, West. Low-coverage single-cell mRNA sequencing reveals cellular heterogeneity and activated signaling pathways in developing cerebral cortex. _Nature Biotechnology_, __32__, 1053-1058 (2014).

Mahata, Zhang, Kolodziejczyk, Proserpio, Haim-Vilmovsky, Taylor, Hebenstreit, Dingler, Moignard, Gottgens, Arlt, McKenzie, Teichmann. Single-Cell RNA Sequencing Reveals T Helper Cells Synthesizing Steroids De Novo to Contribute to Immune Homeostasis. _Cell Reports_, __7__(4): 1130–1142 (2014).

Shalek, Satija, Shuga, Trombetta, Gennert, Lu, Chen, Gertner, Gaublomme, Yosef, Schwartz, Fowler, Weaver, Wang, Wang, Ding, Raychowdhury, Friedman, Hacohen, Park, May, Regev. Single-cell RNA-seq reveals dynamic paracrine control of cellular variation. _Nature_, __510__, 363-369 (2014).

Tasic, Menon, Nguyen, Kim, Jarsky, Yao, Levi, Gray, Sorensen, Dolbeare, Bertagnolli, Goldy, Shapovalova, Parry, Lee, Smith, Bernard, Madisen, Sunkin, Hawrylycz, Koch, Zeng. Adult mouse cortical cell taxonomy revealed by single cell transcriptomics. _Nature Neuroscience_, __19__, 335–346 (2016).

## Disclaimer

The data were downloaded from NCBI's SRA or from EBI's ArrayExpress. See the following disclaimers for copyright / license policies.

http://www.insdc.org/policy.html

http://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html
