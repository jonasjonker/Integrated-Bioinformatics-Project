FROM jonasjonker/mofa2-tutorial:0.0.4

USER root

RUN     apt-get update \ 
    &&  apt-get install -y \
        libhdf5-serial-dev \
    &&  apt-get purge \
    &&  apt-get clean \
    &&  rm -rf /var/lib/apt/lists/* 

USER $NB_USER

RUN     R -e "install.packages('hdf5r')"

RUN     wget --quiet https://bioconductor.org/packages/release/data/annotation/src/contrib/EnsDb.Hsapiens.v86_2.99.0.tar.gz \
    &&  R -e "install.packages('EnsDb.Hsapiens.v86_2.99.0.tar.gz', repos = NULL, type='source')" \
    &&  rm EnsDb.Hsapiens.v86_2.99.0.tar.gz
