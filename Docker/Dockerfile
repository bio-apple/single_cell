FROM ubuntu:latest AS base

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends  libcurl4-openssl-dev\
    build-essential wget curl ca-certificates bzip2 gnupg git libssl-dev libffi-dev libz-dev libbz2-dev && rm -rf /var/lib/apt/lists/* && \
    mkdir -p /ref/ /script/ /raw_data/ /outdir/ && \
    wget -q -O /miniforge.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && \
    bash /miniforge.sh -b -p /opt/conda && rm /miniforge.sh

FROM base AS seurat
RUN /opt/conda/bin/conda create -n seurat --channel conda-forge --channel bioconda --channel defaults r-seurat r-remotes r-signac dnachun::r-cellchat &&  \
    /opt/conda/envs/seurat/bin/Rscript -e "remotes::install_github(c('satijalab/seurat-data', 'satijalab/azimuth', 'satijalab/seurat-wrappers'), quiet = TRUE)"

FROM seurat AS snapy
RUN /opt/conda/bin/conda create -n scanpy --channel conda-forge --channel bioconda --channel defaults scanpy python-igraph leidenalg