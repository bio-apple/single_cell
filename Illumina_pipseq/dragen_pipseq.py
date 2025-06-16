import os
import argparse
import subprocess


def run(pe1,pe2,ref,gtf,outdir,prefix,dragen="/opt/dragen/4.4.4/bin/dragen"):
    pe1=os.path.abspath(pe1)
    pe2=os.path.abspath(pe2)
    ref=os.path.abspath(ref)
    gtf=os.path.abspath(gtf)
    os.makedirs(outdir,exist_ok=True)
    with open(f"{outdir}/{prefix}.fastqlist","w") as f:
        f.write("RGID,RGSM,RGLB,Lane,Read1File,Read2File\n")
        f.write(f"{prefix},{prefix},UnknownLibrary,1,{pe1},{pe2}\n")

    cmd=(f"{dragen} --fastq-list={outdir}/{prefix}.fastqlist "
         f"--ref-dir={ref} --annotation-file={gtf} "
         f"--output-directory={outdir} "
         f"--output-file-prefix={prefix} "
         f"--scrna-enable-pipseq-mode=true "
         f"--fastq-list-sample-id {prefix}")

    print(cmd)
    subprocess.check_call(cmd,shell=True)



if __name__=="__main__":
    parser=argparse.ArgumentParser("This script will run single cell.")
    parser.add_argument("-p1","--pe1",help="R1 fastq file",required=True)
    parser.add_argument("-p2", "--pe2", help="R2 fastq file", required=True)
    parser.add_argument("-o","--outdir",help="directory of output",required=True)
    parser.add_argument("-p","--prefix",help="prefix of output",required=True)
    parser.add_argument("-r","--ref",help="hash build directory",required=True)
    parser.add_argument("-g","--gtf",help="gtf file",required=True)
    parser.add_argument("-d","--dragen",help="path of which dragen",default="/opt/dragen/4.4.4/bin/dragen")
    args=parser.parse_args()
    run(args.pe1, args.pe2, args.ref, args.gtf, args.outdir, args.prefix, args.dragen)