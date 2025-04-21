## Introduction

**PopMAG** is a bioinformatics pipeline that integrates genome-resolved metagenomics with populations genomics workflows. 
The pipeline takes a MAGs samplesheet, a reads files samplesheet and a metadata file as input, and is divided in three main sections:

- MAG Quality control
- Abundance calculation and Variant Calling
- Population genomics metrics and Visualization 



<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

<!--1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))-->
<!--2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))-->

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):-->

First, prepare the samplesheets with your input data that look as follows:

`mag_samplesheet.tsv`:

```csv
sample_id	mag_id	mag_path
SAMPLE_1	CONCOCT_59	MAGs_folder/CONCOCT_59.fa
```

Each row represents a MAG.

`reads_samplesheet.tsv`:

```csv
sample_id	forward	reverse
SAMPLE_1	/forward/reads/path	/reverse/reads/path
```

`Metadata.tsv`:

The Metadata.tsv file must have at least the column sample_id. 

```csv
sample_id	Metadata_1	..	Metadata_n	
SAMPLE_1	0.5	..	0.9
```

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run main.nf \	
   -profile docker \
   --mag_paths mag_samplesheet.csv \
   --reads_paths reads_samplesheet.csv \
   --metadata metadata.csv \ 
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

<!--For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/magenomics/usage) and the [parameter documentation](https://nf-co.re/magenomics/parameters).-->


<!--## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/magenomics/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/magenomics/output).

## Credits

nf-core/magenomics was originally written by Daniel Antonio Sabogal.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->


<!--## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#magenomics` channel](https://nfcore.slack.com/channels/magenomics) (you can join with [this invite](https://nf-co.re/join/slack)).-->

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/magenomics for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->


<!--An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:-->

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
