import scanpy as sc
import os
import datetime
import matplotlib.pyplot as plt
import gzip
import argparse
import base64
import numpy as np
import scipy.sparse as sp
import csv
from io import BytesIO
from dominate import document
from dominate.tags import *

# ------------------- 文件读取 -------------------
def load_adata(prefix):
    adata = sc.read_mtx(f"{prefix}.matrix.mtx.gz").T
    with gzip.open(f"{prefix}.barcodes.tsv.gz", "rt") as f:
        adata.obs_names = [line.strip() for line in f]
    with gzip.open(f"{prefix}.features.tsv.gz", "rt") as f:
        adata.var_names = [line.strip().split('\t')[0] for line in f]
    adata.var_names_make_unique()
    return adata

# ------------------- 图像保存/转换 -------------------
def fig_to_svg_base64(fig):
    buf = BytesIO()
    fig.savefig(buf, format="svg")
    plt.close(fig)
    encoded = base64.b64encode(buf.getvalue()).decode('utf-8')
    return f"data:image/svg+xml;base64,{encoded}"

# ------------------- Knee Plot -------------------
def generate_knee_plot(adata_raw, adata_flt, sample_name):
    raw_sorted = sorted(adata_raw.obs["n_counts"], reverse=True)
    flt_sorted = sorted(adata_flt.obs["n_counts"], reverse=True)
    cutoff_index = len(flt_sorted)
    cutoff_value = raw_sorted[cutoff_index - 1] if cutoff_index <= len(raw_sorted) else raw_sorted[-1]

    fig = plt.figure(figsize=(6, 4))
    plt.plot(range(len(raw_sorted)), raw_sorted, label="Raw", color="gray")
    plt.plot(range(len(flt_sorted)), flt_sorted, label="Filtered", color="blue")
    plt.axvline(x=cutoff_index, linestyle="--", color="red", label=f"Cutoff @ {cutoff_index}")
    plt.text(cutoff_index, cutoff_value, f"{cutoff_index}", color="red")
    plt.xlabel("Barcode Rank")
    plt.ylabel("Total UMI Counts")
    plt.title(f"Knee Plot - {sample_name}")
    plt.yscale("log")
    plt.xscale("log")
    plt.legend()
    plt.tight_layout()
    return fig_to_svg_base64(fig)

# ------------------- Violin Plot -------------------
def generate_violin_plot(adata_flt, sample_name):
    sc.pp.calculate_qc_metrics(adata_flt, inplace=True)
    sc.pl.violin(
        adata_flt,
        ['n_genes_by_counts', 'total_counts'],
        jitter=0.4,
        multi_panel=True,
        show=False
    )
    fig = plt.gcf()
    return fig_to_svg_base64(fig)

# ------------------- Saturation Curve (10x定义) -------------------
def generate_saturation_curve(adata_flt, sample_name, steps=10, seed=42):
    X = adata_flt.X
    if not sp.issparse(X):
        X = sp.coo_matrix(X)
    else:
        X = X.tocoo()

    umi_data = list(zip(X.row, X.col, X.data))
    np.random.seed(seed)
    fractions = np.linspace(0.1, 1.0, steps)
    total_umis_list = []
    unique_umis_list = []

    for frac in fractions:
        sampled = []
        for i, j, count in umi_data:
            k = np.random.binomial(count, frac)
            if k > 0:
                sampled.extend([(i, j)] * k)
        total = len(sampled)
        unique = len(set(sampled))
        total_umis_list.append(total)
        unique_umis_list.append(unique)

    saturation = 1 - (np.array(unique_umis_list) / np.array(total_umis_list))

    fig = plt.figure(figsize=(6, 4))
    plt.plot(total_umis_list, saturation, marker="o", color="purple")
    plt.xlabel("Simulated Sequencing Depth (Total Reads)")
    plt.ylabel("Sequencing Saturation")
    plt.title(f"Saturation Curve - {sample_name}")
    plt.grid(True)
    plt.ylim(0, 1)
    plt.tight_layout()
    return fig_to_svg_base64(fig)

# ------------------- Efficiency Curve (Total UMI vs Median Genes/Cell) -------------------
def generate_efficiency_curve(adata_flt, sample_name, steps=10, seed=42):
    X = adata_flt.X.tocoo()
    umi_data = list(zip(X.row, X.col, X.data))
    np.random.seed(seed)

    fractions = np.linspace(0.1, 1.0, steps)
    sequencing_depth = []
    median_genes = []

    for frac in fractions:
        sampled = []
        for i, j, count in umi_data:
            k = np.random.binomial(count, frac)
            if k > 0:
                sampled.extend([(i, j)] * k)

        rows, cols = zip(*sampled) if sampled else ([], [])
        data = [1] * len(sampled)
        M = sp.coo_matrix((data, (rows, cols)), shape=X.shape).tocsr()

        depth = M.sum()
        genes_per_cell = (M > 0).sum(axis=1).A1
        med_genes = np.median(genes_per_cell) if genes_per_cell.size else 0

        sequencing_depth.append(depth)
        median_genes.append(med_genes)

    fig = plt.figure(figsize=(6, 4))
    plt.plot(sequencing_depth, median_genes, marker="o", color="orange")
    plt.xlabel("Simulated Sequencing Depth (Total UMI)")
    plt.ylabel("Median Genes per Cell")
    plt.title(f"Efficiency Curve - {sample_name}")
    plt.grid(True)
    plt.tight_layout()

    return fig_to_svg_base64(fig)

# ------------------- 保存CSV -------------------
def save_summary_csv(sample_data, output_csv_path):
    header = ["Sample", "Raw_Barcodes", "Filtered_Barcodes", "Genes", "Median_Genes_per_Cell", "Median_UMIs_per_Cell"]
    rows = [[name, d["raw_n"], d["flt_n"], d["genes"], d["med_genes"], d["med_counts"]]
            for name, d in sample_data.items()]

    with open(output_csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)
    print(f"🧾 Summary CSV saved to {output_csv_path}")

# ------------------- HTML 报告 -------------------
def generate_html_report(sample_data, html_path):
    doc = document(title='Single-cell RNA-seq QC Report')

    css = """
    body { font-family: Arial, sans-serif; margin: 40px; background-color: #f9f9f9; }
    h1 { color: #2c3e50; }
    h2 { color: #3a5f8f; margin-top: 60px; }
    table { border-collapse: collapse; width: 100%; margin-bottom: 40px; background: #fff; box-shadow: 0 1px 3px rgba(0,0,0,0.1);}
    th, td { border: 1px solid #ddd; padding: 8px; text-align: center; }
    th { background-color: #e0e5ec; }
    .imgrow { display: flex; flex-wrap: wrap; justify-content: flex-start; gap: 20px; }
    .imgbox {
        background: #fff;
        padding: 10px;
        border-radius: 8px;
        box-shadow: 0 2px 6px rgba(0,0,0,0.1);
        width: 320px;
        text-align: center;
        transition: transform 0.2s ease-in-out;
    }
    .imgbox:hover { transform: scale(1.02); }
    .imgbox img {
        width: 100%;
        height: auto;
        border-radius: 4px;
        cursor: pointer;
        transition: 0.3s;
    }
    .modal {
        display: none;
        position: fixed;
        z-index: 10;
        padding-top: 60px;
        left: 0; top: 0;
        width: 100%; height: 100%;
        background-color: rgba(0,0,0,0.8);
    }
    .modal-content {
        margin: auto;
        display: block;
        max-width: 80%;
        max-height: 80%;
    }
    .modal-content:hover { cursor: zoom-out; }
    """
    js = """
    function enableImageModal() {
        document.querySelectorAll('.imgbox img').forEach(function(img) {
            img.onclick = function() {
                var modal = document.getElementById("imgModal");
                var modalImg = document.getElementById("modalImg");
                modal.style.display = "block";
                modalImg.src = this.src;
            }
        });
        document.getElementById("imgModal").onclick = function() {
            this.style.display = "none";
        }
    }
    window.onload = enableImageModal;
    """

    with doc.head:
        style(css)
        script(js, type="text/javascript")

    with doc:
        h1("Single-cell RNA-seq QC Report")
        p(strong("Generated:"), " ", datetime.datetime.now().strftime('%Y-%m-%d'))

        h2("1. Sample Summary")
        with table():
            with thead():
                with tr():
                    for c in ["Sample", "Raw barcodes", "Filtered barcodes", "Genes", "Median genes/cell", "Median UMIs/cell"]:
                        th(c)
            with tbody():
                for name, info in sample_data.items():
                    with tr():
                        td(name)
                        td(info["raw_n"])
                        td(info["flt_n"])
                        td(info["genes"])
                        td(info["med_genes"])
                        td(info["med_counts"])

        def add_section(title, key):
            h2(title)
            with div(cls="imgrow"):
                for name, info in sample_data.items():
                    with div(cls="imgbox"):
                        h3(name)
                        img(src=info[key])

        add_section("2. Knee Plots", "knee_plot")
        add_section("3. Violin Plots", "violin_plot")
        add_section("4. Saturation Curves", "saturation_plot")
        add_section("5. Efficiency Curves", "efficiency_plot")

        with div(id="imgModal", cls="modal"):
            img(id="modalImg", cls="modal-content")

    os.makedirs(os.path.dirname(html_path), exist_ok=True)
    with open(html_path, "w") as f:
        f.write(doc.render())
    print(f"✅ HTML report saved to {html_path}")

# ------------------- 主函数入口 -------------------
def main():
    parser = argparse.ArgumentParser(description="Generate multi-sample scRNA-seq QC report")
    parser.add_argument("-s", "--samples", nargs='+', required=True, help="List of sample name prefixes")
    parser.add_argument("-b", "--base_dir", nargs='+', required=True, help="List of base directories for each sample")
    parser.add_argument("-o", "--output_dir", default=os.getcwd() + "/scrna_report_output", help="Output directory")
    parser.add_argument("--output_html", default=None, help="Optional: Output HTML filename")
    args = parser.parse_args()

    if len(args.base_dir) != 1 and len(args.samples) != len(args.base_dir):
        raise ValueError("Number of base_dir entries must be 1 or match number of samples.")

    os.makedirs(args.output_dir, exist_ok=True)
    fig_dir = os.path.join(args.output_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    output_html_path = args.output_html or os.path.join(args.output_dir, "multi_sample_scrna_report.html")
    output_csv_path = os.path.join(args.output_dir, "sample_summary.csv")

    sample_data = {}

    for i, sample_name in enumerate(args.samples):
        sample_dir = args.base_dir[0] if len(args.base_dir) == 1 else args.base_dir[i]
        print(f"🔍 Processing sample: {sample_name} in {sample_dir}")

        raw_prefix = os.path.join(sample_dir, f"{sample_name}.scRNA")
        flt_prefix = os.path.join(sample_dir, f"{sample_name}.scRNA.filtered")

        adata_raw = load_adata(raw_prefix)
        adata_flt = load_adata(flt_prefix)

        adata_raw.obs["n_counts"] = adata_raw.X.sum(axis=1).A1
        adata_raw.obs["n_genes"] = (adata_raw.X > 0).sum(axis=1).A1
        adata_flt.obs["n_counts"] = adata_flt.X.sum(axis=1).A1
        adata_flt.obs["n_genes"] = (adata_flt.X > 0).sum(axis=1).A1

        knee_base64 = generate_knee_plot(adata_raw, adata_flt, sample_name)
        violin_base64 = generate_violin_plot(adata_flt, sample_name)
        saturation_base64 = generate_saturation_curve(adata_flt, sample_name)
        efficiency_base64 = generate_efficiency_curve(adata_flt, sample_name)

        sample_data[sample_name] = {
            'raw_n': adata_raw.n_obs,
            'flt_n': adata_flt.n_obs,
            'genes': adata_raw.n_vars,
            'med_genes': int(np.median(adata_flt.obs['n_genes'])),
            'med_counts': int(np.median(adata_flt.obs['n_counts'])),
            'knee_plot': knee_base64,
            'violin_plot': violin_base64,
            'saturation_plot': saturation_base64,
            'efficiency_plot': efficiency_base64
        }

    generate_html_report(sample_data, output_html_path)
    save_summary_csv(sample_data, output_csv_path)

if __name__ == "__main__":
    main()
