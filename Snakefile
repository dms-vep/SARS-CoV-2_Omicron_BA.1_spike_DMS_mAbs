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
        clade_desc=lambda wc: (
            f"of all clades" if wc.clade == "all" else f"of {wc.clade}"
        ),
    output:
        nb="results/compare_muteffects/compare_muteffects_{clade}.ipynb",
        csv="results/compare_muteffects/{clade}_natural_enrichment_vs_dms.csv",
        html="results/compare_muteffects/{clade}_natural_enrichment_vs_dms.html",
    log:
        os.path.join(config["logdir"], "compare_muteffects_{clade}.txt"),
    conda:
        "dms-vep-pipeline/environment.yml"
    shell:
        """
        papermill -p clade {wildcards.clade} {input.nb} {output.nb} &> {log}
        python {input.pyscript} \
            --chart {output.html} \
            --output {output.html} \
            --markdown {input.md} \
            --site {params.site} \
            --title \
            "DMS vs natural evolution {params.clade_desc} for {config[github_repo]}" \
            --description \
            "Correlations among DMS datasets and natural evolution {params.clade_desc}" \
            &>> {log}
        """


# Add any extra data/results files for docs with name: file
extra_data_files = {
    "sequential to reference site numbering": config["site_numbering_map"],
}
extra_data_files.update(
    {
        f"{clade} clade natural evolution enrichment vs DMS": (
            f"results/compare_muteffects/{clade}_natural_enrichment_vs_dms.csv"
        )
        for clade in config["actual_vs_expected_clades"]
    }
)

# Add any extra HTML documents to display here with name: file
extra_html_docs = {
    f"{clade} clade natural enrichment versus DMS": (
        f"results/compare_muteffects/{clade}_natural_enrichment_vs_dms.html"
    )
    for clade in config["actual_vs_expected_clades"]
}

# If you add rules with "nb" output that have wildcards, specify the rule name
# and subindex titles for the wildcards as in "docs.smk" for `nb_rule_wildcards`
# and `subindex_titles`
extra_nb_rule_wildcards = {
    "compare_muteffects": {"clade": config["actual_vs_expected_clades"]}
}
extra_subindex_titles = {"compare_muteffects": "natural evolution enrichment vs DMS"}

# include `dms-vep-pipeline` docs building Snakemake file
include: os.path.join(config["pipeline_path"], "docs.smk")
