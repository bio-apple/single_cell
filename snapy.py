import scanpy as sc
import os
import datetime
import matplotlib.pyplot as plt
import gzip
import argparse
import base64
from io import BytesIO

def load_adata(prefix):
    adata = sc.read_mtx(f"{prefix}.matrix.mtx.gz").T
    with gzip.open(f"{prefix}.barcodes.tsv.gz", "rt") as f:
        adata.obs_names = [line.strip() for line in f]
    with gzip.open(f"{prefix}.features.tsv.gz", "rt") as f:
        adata.var_names = [line.strip().split('\t')[0] for line in f]
    adata.var_names_make_unique()
    return adata

def fig_to_svg_base64(fig):
    buf = BytesIO()
    fig.savefig(buf, format="svg")
    plt.close(fig)
    encoded = base64.b64encode(buf.getvalue()).decode('utf-8')
    return f"data:image/svg+xml;base64,{encoded}"

def save_svg_png(fig, base_path):
    svg_path = base_path + ".svg"
    png_path = base_path + ".png"
    fig.savefig(svg_path)
    fig.savefig(png_path, dpi=300)
    encoded_svg = fig_to_svg_base64(fig)
    return encoded_svg


def generate_knee_plot(adata_raw, adata_flt, sample_name, fig_dir):
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
    return fig_to_svg_base64(fig)  # 返回Base64字符串


def generate_violin_plot(adata_flt, sample_name, fig_dir):
    sc.pp.calculate_qc_metrics(adata_flt, inplace=True)
    base_name = f"{sample_name}_violin"

    sc.pl.violin(
        adata_flt,
        ['n_genes_by_counts', 'total_counts'],
        jitter=0.4,
        multi_panel=True,
        show=False
    )

    fig = plt.gcf()
    return fig_to_svg_base64(fig)  # 返回Base64字符串


def generate_html_report(sample_data, html_path):
    summary_rows = []
    for name, info in sample_data.items():
        summary_rows.append(
            f"<tr><td>{name}</td><td>{info['raw_n']}</td><td>{info['flt_n']}</td><td>{info['genes']}</td>"
            f"<td>{info['med_genes']}</td><td>{info['med_counts']}</td></tr>")

    html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Multi-sample scRNA Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        h1 {{ color: #2c3e50; }}
        table {{ border-collapse: collapse; width: 100%; margin-bottom: 40px; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: center; }}
        th {{ background-color: #f2f2f2; }}
        .imgrow {{ display: flex; flex-direction: row; gap: 20px; flex-wrap: wrap; margin-bottom: 40px; }}
        .imgbox {{ border: 1px solid #ccc; padding: 10px; text-align: center; }}
        .imgbox img {{ max-width: 400px; height: auto; }}
    </style>
</head>
<body>
    <h1>Single-cell RNA-seq Report</h1>
    <p><strong>Date:</strong> {datetime.datetime.now().strftime('%Y-%m-%d')}</p>

    <h2>1. Sample Summary Table</h2>
    <table>
        <tr>
            <th>Sample</th>
            <th>Total Raw Barcodes</th>
            <th>Total Filtered Barcodes</th>
            <th>Total Genes</th>
            <th>Median Genes/Cell (Filtered)</th>
            <th>Median UMIs/Cell (Filtered)</th>
        </tr>
        {''.join(summary_rows)}
    </table>

    <h2>2. Knee Plots</h2>
    <div class="imgrow">
        {''.join([f'<div class="imgbox"><h3>{name}</h3><img src="{info["knee_plot"]}"/></div>' for name, info in sample_data.items()])}
    </div>

    <h2>3. Violin Plots</h2>
    <div class="imgrow">
        {''.join([f'<div class="imgbox"><h3>{name}</h3><img src="{info["violin_plot"]}"/></div>' for name, info in sample_data.items()])}
    </div>
</body>
</html>
"""

    os.makedirs(os.path.dirname(html_path), exist_ok=True)
    with open(html_path, "w") as f:
        f.write(html_content)
    print(f"✅ Report saved to {html_path}")

def main():
    parser = argparse.ArgumentParser(description="Generate multi-sample scRNA-seq QC report")
    parser.add_argument("-s",'--samples', nargs='+', required=True, help="List of sample name prefixes")
    parser.add_argument("-b",'--base_dir', nargs='+', required=True, help="List of base directories corresponding to each sample")
    parser.add_argument('-o','--output_dir', help=f"Output folder for HTML and figures",default=os.getcwd()+"/scrna_report_output")
    parser.add_argument('-html','--output_html', default=None, help="Optional: Output HTML filename (default: multi_sample_scrna_report.html in output_dir)")
    args = parser.parse_args()

    if len(args.base_dir) != 1 and len(args.samples) != len(args.base_dir):
        raise ValueError("The number of base_dir entries must be 1 or match the number of samples.")

    os.makedirs(args.output_dir, exist_ok=True)
    fig_dir = os.path.join(args.output_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    output_html_path = args.output_html or os.path.join(args.output_dir, "multi_sample_scrna_report.html")

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

        knee_base64 = generate_knee_plot(adata_raw, adata_flt, sample_name, fig_dir)
        violin_base64 = generate_violin_plot(adata_flt, sample_name, fig_dir)

        sample_data[sample_name] = {
            'raw_n': adata_raw.n_obs,
            'flt_n': adata_flt.n_obs,
            'genes': adata_raw.n_vars,
            'med_genes': int(adata_flt.obs['n_genes'].median()),
            'med_counts': int(adata_flt.obs['n_counts'].median()),
            'knee_plot': knee_base64,
            'violin_plot': violin_base64,
        }

    generate_html_report(sample_data, output_html_path)


if __name__ == "__main__":
    main()
