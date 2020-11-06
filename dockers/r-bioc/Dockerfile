FROM rocker/r-base:4.0.3

LABEL maintainer="jonasjonker@gmail.com"

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
        curl \
        git \
        libcairo2-dev \
        libcurl4-openssl-dev \
        libfreetype6-dev \
        libfribidi-dev \
        libgit2-dev \
        libgsl-dev \
        libgtk2.0-dev \
        libharfbuzz-dev \
        libhdf5-serial-dev \
        libjpeg-dev \
        libpng-dev \
        libssl-dev \
        libtiff5-dev \
        libxml2-dev \
        libxt-dev \
        openssl \
        pkg-config \
        python3-dev \
        python3-pip \
        python3-venv \
        zlib1g-dev \
    &&  apt-get purge \
    &&  apt-get clean \
    &&  rm -rf /var/lib/apt/lists/* 

RUN     mkdir -p ${VENV_DIR} \
    &&  chown -R ${NB_USER} ${VENV_DIR}

USER    ${NB_USER}

RUN     python3 -m venv ${VENV_DIR} \
    &&  pip3 install pip==20.2.3 \ 
    &&  pip3 install --no-cache-dir \
        jupyter-rsession-proxy \
    &&  pip3 install mofapy2==0.5.6

RUN     R -e "install.packages('Cairo')"

RUN     R -e "install.packages(c('devtools','curl'))"  \
    &&  R -e "install.packages(c('tidyverse','docopt', 'testthat'))"

RUN     echo "export PATH=${PATH}" >> ${HOME}/.profile

RUN     R -e "devtools::install_github('RBigData/pbdZMQ')" \
    &&  R -e "devtools::install_github('IRkernel/IRkernel')" \
    &&  R -e "IRkernel::installspec(prefix='${VENV_DIR}')"

RUN     R -e "install.packages('BiocManager')" 

CMD jupyter notebook --ip 0.0.0.0
