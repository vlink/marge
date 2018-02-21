# Copyright Verena M. Link <vlink@ucsd.edu>
# 
# This dockerfile was generated to install MARGE
#
# MARGE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MARGE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

#Get ubuntu
FROM ubuntu:14.04
RUN apt-get -y -m update && apt-get -y -m install apt-utils g++ gcc make zip unzip wget vim 

RUN echo 'export PATH=$PATH:/opt/conda/bin' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-4.3.14-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

RUN apt-get install -y curl grep sed dpkg && \
    TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb && \
    apt-get clean

ENV PATH $PATH:/opt/conda/bin

##Install conda, bioconda and some tools necessary for HOMER
RUN conda config --add channels conda-forge
RUN conda config --add channels r
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda install r-essentials bioconductor-deseq2 bioconductor-edger r-seqinr rstudio

##Install HOMER
RUN wget http://homer.ucsd.edu/homer/configureHomer.pl
RUN mkdir homer && mv configureHomer.pl homer
RUN perl homer/configureHomer.pl -install
ENV PATH $PATH://homer//bin/
#Add path to bashrc

##Install MARGE
RUN wget http://homer.ucsd.edu/MARGE/MARGE_v1.0.tar.gz
RUN tar fxvz MARGE_v1.0.tar.gz
RUN chown -R root:root /MARGE_v1.0*
RUN curl -L http://cpanmin.us | perl - --self-upgrade
RUN cpanm Set::IntervalTree
RUN cpanm Statistics::Basic
RUN cpanm Sys::CPU
WORKDIR /MARGE_v1.0
RUN perl configure.pl
RUN perl install.pl
ENV PATH $PATH:/MARGE_v1.0/bin/
RUN chmod -R a+x /MARGE_v1.0/bin/
WORKDIR /
