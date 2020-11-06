FROM jonasjonker/r-bioc:0.0.1

RUN     R -e "devtools::install_github('GreenleafLab/ArchR', ref='master', repos = BiocManager::repositories())" \
    &&  R -e "ArchR::installExtraPackages()"

RUN     wget --quiet http://bioconductor.org/packages/release/data/annotation/src/contrib/BSgenome.Hsapiens.UCSC.hg38_1.4.3.tar.gz \
    &&  R -e "install.packages('BSgenome.Hsapiens.UCSC.hg38_1.4.3.tar.gz', repos = NULL, type='source')" \
    &&  rm BSgenome.Hsapiens.UCSC.hg38_1.4.3.tar.gz

CMD jupyter notebook --ip 0.0.0.0
