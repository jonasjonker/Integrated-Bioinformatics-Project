FROM jonasjonker/mofa2-notebook:0.0.2 AS baselayer

USER    root

RUN     apt-get update \ 
    &&  apt-get install -y \
        pkg-config \
        libgsl-dev \
        zlib1g-dev \
        libgtk2.0-dev \
    &&  apt-get purge \
    &&  apt-get clean \
    &&  rm -rf /var/lib/apt/lists/* 

USER    $NB_USER

FROM baselayer AS layer-2

RUN     R -e "install.packages('BiocManager')"

FROM layer-2 AS layer-3

RUN     R -e "BiocManager::install('motifmatchr')"

FROM layer-3 AS layer-4

RUN     R -e "BiocManager::install('Signac')"

FROM layer-4 AS layer-5

RUN     R -e "BiocManager::install(c('msigdbr', 'TFBSTools', 'JASPAR2020'))"

FROM layer-5 AS layer-6

RUN     wget --quiet http://bioconductor.org/packages/release/data/annotation/src/contrib/BSgenome.Hsapiens.UCSC.hg38_1.4.3.tar.gz \
    &&  R -e "install.packages('BSgenome.Hsapiens.UCSC.hg38_1.4.3.tar.gz', repos = NULL, type='source')" \
    &&  rm BSgenome.Hsapiens.UCSC.hg38_1.4.3.tar.gz

CMD jupyter notebook --ip 0.0.0.0
