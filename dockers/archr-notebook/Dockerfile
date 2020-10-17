FROM rocker/r-base:4.0.3 AS base-layer

ENV     NB_USER docker
ENV     NB_UID 1000
ENV     NB_GID 1000
ENV     VENV_DIR /srv/venv
ENV     HOME /home/${NB_USER}

ENV LD_LIBRARY_PATH /usr/local/lib/R/lib
ENV PATH ${VENV_DIR}/bin:$PATH

WORKDIR ${HOME}

RUN     apt-get update \ 
    &&  apt-get install -y \
        libcairo2-dev \
        libxt-dev \
        python3-venv \
        python3-dev \
        python3-pip \
        libcurl4-openssl-dev \
        curl \
        libssl-dev \
        libgit2-dev \
        openssl \
        libxml2-dev \
        git \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
    &&  apt-get purge \
    &&  apt-get clean \
    &&  rm -rf /var/lib/apt/lists/* 

RUN     mkdir -p ${VENV_DIR} \
    &&  chown -R ${NB_USER} ${VENV_DIR}

USER    ${NB_USER}

RUN     python3 -m venv ${VENV_DIR} \
    &&  pip3 install pip==20.2.3 \ 
    &&  pip3 install --no-cache-dir \
        jupyter-rsession-proxy

RUN     R -e "install.packages('Cairo')"

RUN     R -e "install.packages(c('devtools','curl'))" 

RUN     echo "export PATH=${PATH}" >> ${HOME}/.profile

RUN     R -e "devtools::install_github('RBigData/pbdZMQ')" \
    &&  R -e "devtools::install_github('IRkernel/IRkernel')" \
    &&  R -e "IRkernel::installspec(prefix='${VENV_DIR}')"

FROM base-layer AS bioc-layer

RUN     R -e "install.packages('BiocManager')"

FROM bioc-layer AS archer-layer

USER root

RUN     apt-get update \ 
    &&  apt-get install -y \
        libgsl-dev \
    &&  apt-get purge \
    &&  apt-get clean \
    &&  rm -rf /var/lib/apt/lists/* 

USER $NB_USER

RUN     R -e "devtools::install_github('GreenleafLab/ArchR', ref='master', repos = BiocManager::repositories())"

RUN     R -e "ArchR::installExtraPackages()"

FROM archer-layer AS reference-layer

RUN     wget --quiet http://bioconductor.org/packages/release/data/annotation/src/contrib/BSgenome.Hsapiens.UCSC.hg38_1.4.3.tar.gz \
    &&  R -e "install.packages('BSgenome.Hsapiens.UCSC.hg38_1.4.3.tar.gz', repos = NULL, type='source')" \
    &&  rm BSgenome.Hsapiens.UCSC.hg38_1.4.3.tar.gz

CMD jupyter notebook --ip 0.0.0.0
