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

