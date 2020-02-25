FROM continuumio/anaconda:5.3.0

# Install IGoR version 1.3.0
RUN apt-get update -y && \
  apt-get install -y build-essential unzip python3-pip && \
  wget https://github.com/qmarcou/IGoR/releases/download/1.3.0/igor_1-3-0.zip && \
  unzip igor_1-3-0.zip
WORKDIR /igor_1-3-0
RUN ./configure && \
  make && \
  make install 
RUN pip3 install ./pygor
RUN conda update -y conda
RUN conda create -n pygor python=3 pandas biopython matplotlib numpy scipy
ENV PYGOR_PATH="/igor_1-3-0/pygor"
WORKDIR ..

# Install Change-O
RUN apt-get install idle3 -y && \
  pip3 install --no-cache-dir --upgrade pip

# Delete cache and tmp files
RUN apt-get clean && \
  apt-get autoclean && \
  rm -rf /var/cache/* && \
  rm -rf /tmp/* && \
  rm -rf /var/tmp/* && \
  rm -rf /var/lib/apt/lists/* 

RUN pip3 install hg+https://bitbucket.org/kleinstein/changeo#default --user
ENV CHANGEO_DIR="~/.local/bin"

# Install IgBlast version 1.13.0
RUN wget ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/1.13.0/ncbi-igblast-1.13.0-x64-linux.tar.gz 
RUN tar -zxf ncbi-igblast-1.13.0-x64-linux.tar.gz
ENV IGBLAST_DIR="/ncbi-igblast-1.13.0"
ENV PATH="${IGBLAST_DIR}/bin:${PATH}"

# This builds the IMGT database for IgBlast + Change-O
RUN wget https://bitbucket.org/kleinstein/immcantation/raw/3.0.0/scripts/fetch_igblastdb.sh
RUN mv fetch_igblastdb.sh /usr/local/bin

RUN wget https://bitbucket.org/kleinstein/immcantation/raw/3.0.0/scripts/fetch_imgtdb.sh
RUN mv fetch_imgtdb.sh /usr/local/bin

RUN wget https://bitbucket.org/kleinstein/immcantation/raw/3.0.0/scripts/clean_imgtdb.py
RUN mv clean_imgtdb.py /usr/local/bin
RUN chmod +x /usr/local/bin/clean_imgtdb.py

RUN wget https://bitbucket.org/kleinstein/immcantation/raw/3.0.0/scripts/imgt2igblast.sh
RUN mv imgt2igblast.sh /usr/local/bin

RUN bash fetch_igblastdb.sh -o ~/share/igblast
RUN bash fetch_imgtdb.sh -o ~/share/germlines/imgt
RUN bash imgt2igblast.sh -i ~/share/germlines/imgt -o ~/share/igblast
RUN cp -r ~/share/igblast/* ${IGBLAST_DIR}
RUN cp -r ~/share/germlines ${IGBLAST_DIR}

#  partis
RUN git clone https://github.com/psathyrella/partis.git

RUN apt-get update && apt-get install -y \
  build-essential \
    libboost-dev \
    libgsl-dev \
    libncurses-dev \
    libxt6 \
    libyaml-cpp-dev \
    libyaml-dev \
    libz-dev

RUN conda install -y -cbioconda -cbiocore -cr \
  python=2.7 \
  biopython \
  pandas \
  psutil \
  pysam \
  scons \
  seaborn \
  zlib \
  pyyaml \
  scikit-learn \
  pysam \
  mafft \
  r-essentials \
  r-devtools \
  r-roxygen2

RUN pip install \
  colored-traceback \
  dendropy==4.0.0
COPY . /partis
WORKDIR /partis
RUN ./bin/build.sh
ENV PARTIS_PATH="/partis/bin/partis"
WORKDIR ..

RUN unset R_LIBS_SITE
RUN Rscript -e 'install.packages(c("TreeSim", "TreeSimGM"), repos="http://cran.rstudio.com")'

# Download sumrep
RUN git clone --single-branch --branch master https://github.com/matsengrp/sumrep.git /sumrep
ENV SUMREP_PATH="/sumrep"

# Download partis annotations for MDS example
RUN mkdir -p /sumrep/Examples/flu && \
  wget --no-check-certificate https://zenodo.org/record/3385364/files/flu_rds.tar && \
  tar -C /sumrep/Examples/flu -xvf flu_rds.tar && \
  rm flu_rds.tar

# Install sumrep
ARG R_DEPS="c('Rcpp', 'devtools', 'roxygen2', 'testthat', 'rmarkdown', 'knitr', 'optparse', 'BiocManager')"
RUN Rscript -e "install.packages(${R_DEPS}, repos='http://cran.rstudio.com')" && \
  Rscript -e "setwd('/sumrep'); library(devtools); install_deps(dependencies=TRUE, upgrade=TRUE)" && \
  Rscript -e "setwd('/sumrep'); library(devtools); document(); install(build_vignettes=TRUE)"
