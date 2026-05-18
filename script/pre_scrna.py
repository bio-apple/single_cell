import os,sys,re
import subprocess
import argparse

docker="scrna:latest"

def cellbender(h5_file,outdir,prefix,expected_cells = 5000):
    #CellBender:https://github.com/broadinstitute/CellBender
    #CellBender is a software package for eliminating technical artifacts from high-throughput single-cell RNA sequencing (scRNA-seq) data.
    #Fleming S J, Chaffin M D, Arduini A, et al. Unsupervised removal of systematic background noise from droplet-based single-cell experiments using CellBender[J]. Nature methods, 2023, 20(9): 1323-1335.
    h5_file=os.path.abspath(h5_file)
    outdir=os.path.abspath(outdir)
    if not os.path.exists(outdir):
        subprocess.check_call(f'mkdir -p {outdir}',shell=True)
    cmd=(f'docker -v {os.path.dirname(h5_file)}:/raw_data/ -v {outdir}:/outdir {docker} sh -c\'export PATH=/opt/conda/envs/python3.7/bin/:$PATH && '
         f'cellbender remove-background --expected-cells {expected_cells} --input /raw_data/{h5_file.split("/")[-1]} --output /outdir/{prefix}.cellbender_filtered.h5\'')
    subprocess.check_call(cmd,shell=True)

def scDblFinder(input_h5,prefix,outdir):
    #scDblFinder:https://github.com/plger/scDblFinder
    #The scDblFinder package gathers various methods for the detection and handling of doublets/multiplets in single-cell sequencing data (i.e. multiple cells captured within the same droplet or reaction volume), including the novel scDblFinder method.
    #Germain P L, Lun A, Meixide C G, et al. Doublet identification in single-cell sequencing data using scDblFinder[J]. f1000research, 2022, 10: 979.
    input_h5=os.path.abspath(input_h5)
    r_script = f"""
        library(scDblFinder)

        # 读取10X h5
        sce <- read10xCounts("/raw_data/{input_h5.split("/")[-1]}")

        # 运行 scDblFinder
        sce <- scDblFinder(sce)

        # 保留 singlet
        sce_singlet <- sce[, sce$scDblFinder.class == "singlet"]

        # 保存为 RDS
        saveRDS(sce_singlet, "/outdir/{prefix}.filtered_singlet.rds")
        """
    with open(f'{outdir}/run_scdblfinder.R', "w") as f:
        f.write(r_script)


    # 使用 Docker 运行 R 脚本
    cmd = [
        "docker", "run",
        "-v", f"{outdir}:/outdir",
        "-v", f"{os.path.dirname(input_h5)}:/raw_data/",
        docker,
        "/opt/conda/envs/r-base4.4.3/bin/Rscript",
        "/outdir/run_scdblfinder.R"
    ]
    print(f"Running scDblFinder in Docker: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    parser=argparse.ArgumentParser("")
    parser.add_argument("-h5","--h5_file",help="h5 file",required=True)
    parser.add_argument("-e","--expect",help="expect cell number",default=5000,required=True)
    parser.add_argument("-o","--outdir",help="output directory",default=os.getcwd())
    parser.add_argument("-p","--prefix",help="prefix of output",required=True)
    args=parser.parse_args()
    cellbender(args.h5_file,args.outdir,args.prefix,args.expected)
    scDblFinder(args.outdir+"/"+args.prefix+".cellbender_filtered.h5",args.outdir,args.prefix)