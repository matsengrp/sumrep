FROM continuumio/anaconda

# Install Change-O
RUN apt-get update -y && \
  apt-get install python3-pip idle3 -y && \
  pip3 install --no-cache-dir --upgrade pip && \
  \
  # delete cache and tmp files
  apt-get clean && \
  apt-get autoclean && \
  rm -rf /var/cache/* && \
  rm -rf /tmp/* && \
  rm -rf /var/tmp/* && \
  rm -rf /var/lib/apt/lists/* 

RUN pip3 install hg+https://bitbucket.org/kleinstein/changeo#default --user
ENV CHANGEO_DIR="~/.local/bin"

# Install IgBlast
RUN wget ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/ncbi-igblast-1.13.0-x64-linux.tar.gz \
  && tar -zxf ncbi-igblast-1.13.0-x64-linux.tar.gz
ENV IGBLAST_DIR="/ncbi-igblast-1.13.0"
ENV PATH="${IGBLAST_DIR}/bin:${PATH}"

# This builds the IMGT database for IgBlast + Change-O
RUN wget https://bitbucket.org/kleinstein/immcantation/raw/41165c9b9cd10ddf8f27ec80fae799491d6db2db/scripts/fetch_igblastdb.sh
RUN mv fetch_igblastdb.sh /usr/local/bin

RUN wget https://bitbucket.org/kleinstein/immcantation/raw/41165c9b9cd10ddf8f27ec80fae799491d6db2db/scripts/fetch_imgtdb.sh
RUN mv fetch_imgtdb.sh /usr/local/bin

RUN wget https://bitbucket.org/kleinstein/immcantation/raw/41165c9b9cd10ddf8f27ec80fae799491d6db2db/scripts/clean_imgtdb.py
RUN mv clean_imgtdb.py /usr/local/bin
RUN chmod +x /usr/local/bin/clean_imgtdb.py

RUN wget https://bitbucket.org/kleinstein/immcantation/raw/41165c9b9cd10ddf8f27ec80fae799491d6db2db/scripts/imgt2igblast.sh
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

RUN conda update -y conda
RUN conda install -y -cbioconda -cbiocore -cr \
  python=2.7 \
  biopython \
  pandas \
  psutil \
  pysam \
  scons \
  seaborn \
  zlib \
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

RUN unset R_LIBS_SITE
RUN R --vanilla --slave -e 'install.packages(c("TreeSim", "TreeSimGM"), repos="http://cran.rstudio.com/")'

# Next let's install sumrep dependencies
RUN R --vanilla --slave -e \
  'install.packages(c("alakazam", "ape", "CollessLike", "data.table", "dplyr", "entropy", "HDMD", "jsonlite", "magrittr", "Peptides", "RecordLinkage", "shazam", "seqinr", "stringdist", "stringr", "testthat", "textmineR", "yaml"), repos = "http://cran.us.r-project.org")' && \
  R --vanilla --slave -e 'source("https://bioconductor.org/biocLite.R"); biocLite("Biostrings")'

WORKDIR ..
RUN git clone https://github.com/matsengrp/sumrep.git
COPY . /sumrep
