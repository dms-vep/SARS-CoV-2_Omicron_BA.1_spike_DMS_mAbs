"""Top-level ``snakemake`` file that runs analysis."""


import os


configfile: "config.yaml"


# include `dms-vep-pipeline` pipeline Snakemake file
include: os.path.join(config["pipeline_path"], "pipeline.smk")


rule all:
    input:
        variant_count_files,
        rules.check_adequate_variant_counts.output.passed,
        antibody_escape_files,
        (
            [config["muteffects_observed"], config["muteffects_latent"]]
            if len(func_selections)
            else []
        ),
        config["docs"],


# Arbitrary other rules should be added here
rule site_numbering_map:
    """Map sequential numbering of protein in experiments to standard reference."""
    input:
        prot=config["gene_sequence_protein"],
        reference_site_regions=config["reference_site_regions"],
    output:
        reference="results/site_numbering/numbering_reference.fa",
        alignment="results/site_numbering/alignment.fa",
        to_align="results/site_numbering/to_align.fa",
        site_numbering_map=config["site_numbering_map"],
    params:
        numbering_reference_accession=config["numbering_reference_accession"],
    log:
        os.path.join(config["logdir"], "site_numbering_map.txt"),
    conda:
        "dms-vep-pipeline/environment.yml"
    script:
        "scripts/site_numbering_map.py"


rule validation_IC50s:
    """Get ``polyclonal`` predicted IC50s for validated mutations."""
    input:
        config["validation_ic50s"],
        [
            os.path.join(config["escape_dir"], f"{antibody}.pickle")
            for antibody in pd.read_csv(config["validation_ic50s"])["antibody"].unique()
        ],
        nb="notebooks/validation_IC50s.ipynb",
    output:
        nb="results/notebooks/validation_IC50s.ipynb",
    log:
        os.path.join(config["logdir"], "validation_IC50s.txt"),
    conda:
        "dms-vep-pipeline/environment.yml"
    shell:
        "papermill {input.nb} {output.nb} &> {log}"


rule compare_muteffects:
    """Compare measured mutational effects to natural and Wu data."""
    input:
        config["actual_vs_expected_mut_counts"],
        config["muteffects_observed"],
        nb="notebooks/compare_muteffects.ipynb",
        pyscript=os.path.join(config["pipeline_path"], "scripts/format_altair_html.py"),
        md="plot_legends/compare_muteffects.md",
    params:
        config["prefusion_Tan_excel"],
        config["ntd_Ouyang"],
        config["rbd_dms_Starr"],
        config["rbd_dms_Starr_target"],
        config["actual_vs_expected_clades"],
        config["actual_vs_expected_pseudocount"],
        config["actual_vs_expected_min_expected"],
        config["muteffects_plot_kwargs"]["addtl_slider_stats"]["times_seen"],
        site=lambda _, output: os.path.join(
            github_pages_url,
            os.path.basename(output.html),
        ),
    output:
        nb="results/compare_muteffects/compare_muteffects.ipynb",
        csv=config["natural_enrichment_vs_dms"],
        html=config["natural_enrichment_vs_dms_plot"],
    log:
        os.path.join(config["logdir"], "compare_muteffects.txt"),
    conda:
        "dms-vep-pipeline/environment.yml"
    shell:
        """
        papermill {input.nb} {output.nb} &> {log}
        python {input.pyscript} \
            --chart {output.html} \
            --output {output.html} \
            --markdown {input.md} \
            --site {params.site} \
            --title "DMS vs natural evolution for {config[github_repo]}" \
            --description "Correlations among DMS datasets and natural evolution" \
            &>> {log}
        """


# Add any extra data/results files for docs with name: file
extra_data_files = {
    "sequential to reference site numbering": config["site_numbering_map"],
    "natural evolution enrichment vs DMS": config["natural_enrichment_vs_dms"],
}

# Add any extra HTML documents to display here with name: file
extra_html_docs = {
    "natural enrichment versus DMS": config["natural_enrichment_vs_dms_plot"]
}

# include `dms-vep-pipeline` docs building Snakemake file
include: os.path.join(config["pipeline_path"], "docs.smk")
