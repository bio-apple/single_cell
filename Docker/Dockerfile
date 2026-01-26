FROM ubuntu:latest
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends  libcurl4-openssl-dev\
    build-essential wget curl ca-certificates bzip2 gnupg git libssl-dev libffi-dev libz-dev libbz2-dev && rm -rf /var/lib/apt/lists/* && \
    mkdir -p /ref/ /script/ /raw_data/ /outdir/ && \
    wget -q -O /miniforge.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && \
    bash /miniforge.sh -b -p /opt/conda && rm /miniforge.sh &&  \
    /opt/conda/bin/conda create -n r-base4.4.3 --channel conda-forge --channel bioconda --channel defaults r-base==4.4.3 r-seurat==5.3.0 r-remotes r-signac dnachun::r-cellchat &&  \
    /opt/conda/envs/r-base4.4.3/bin/Rscript -e "remotes::install_github(c('satijalab/seurat-data', 'satijalab/azimuth', 'satijalab/seurat-wrappers'), quiet = TRUE)" && \
    /opt/conda/bin/conda create -n r-base4.0.3 --channel conda-forge --channel bioconda --channel defaults --channel bioturing r-base==4.0.3 r-nichenetr  \
    r-devtools r-essentials r-biocmanager bioconductor-clusterprofiler bioconductor-scdblfinder bioconductor-singler r-monocle3 bioconductor-deseq2 bioconductor-aucell celltypist && \
    /opt/conda/bin/conda create -n python3.7 --channel conda-forge --channel bioconda --channel defaults cellbender==0.3.2 && \
    /opt/conda/bin/conda clean -a -y

