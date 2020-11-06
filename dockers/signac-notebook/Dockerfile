FROM jonasjonker/r-bioc:0.0.1

RUN     R -e "install.packages('hdf5r')" \
    &&  R -e "BiocManager::install('motifmatchr')" \
    &&  R -e "BiocManager::install('Signac')" \
    &&  R -e "BiocManager::install(c('msigdbr', 'TFBSTools'))"

RUN     wget --quiet https://bioconductor.org/packages/release/data/annotation/src/contrib/EnsDb.Hsapiens.v86_2.99.0.tar.gz \
    &&  R -e "install.packages('EnsDb.Hsapiens.v86_2.99.0.tar.gz', repos = NULL, type='source')" \
    &&  rm EnsDb.Hsapiens.v86_2.99.0.tar.gz
