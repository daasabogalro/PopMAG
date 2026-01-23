## Introduction

**PopMAG** is a pipeline that integrates genome-resolved metagenomics data with population genomics tools to analyze metagenome-assembled genomes (MAGs) and their population-level variations. The pipeline processes MAGs alongside paired-end sequencing short reads to perform quality assessment, abundance profiling, variant calling, and population genomics analyses, ending in an interactive visualization dashboard built with shiny.

The pipeline is organized into five main phases:

- MAG Quality control and preprocessing.
- Microbial community profiling.
- Abundance calculation and variant calling.
- Population genomics and functional analysis.
- Visualization and reporting.

**PopMAG** can aid to understand both the functional potential and population dynamics of metagenome-assembled genomes, particularly in the context of comparative genomics and temporal or spatial studies.

Complete documentation for the pipeline can be found in [daasabogalro.github.io/popmag_docs/](https://daasabogalro.github.io/popmag_docs/)

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
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. <!--Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.-->

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

`Metadata.csv`:

The Metadata.csv file must have at least the column sample_id.

```csv
sample_id,Metadata_1,Metadata_n
SAMPLE_1,0.5,0.9
```

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run main.nf \
   -profile docker \
   --mag_paths mag_samplesheet.tsv \
   --reads_paths reads_samplesheet.tsv \
   --metadata metadata.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/magenomics/usage) and the [parameter documentation](https://daasabogalro.github.io/popmag_docs/config).

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

<!--An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
-->

> **CheckM2: A rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. (2023).**
>
> Chklovski, A., Parks, D. H., Woodcroft, B. J., & Tyson, G. W.
>
> _Nature Methods_, 20(8), 1203–1212. https://doi.org/10.1038/s41592-023-01940-w

> **dRep: A tool for fast and accurate genomic comparisons that enables improved genome recovery from metagenomes through de-replication. (2017).**
>
> Olm, M. R., Brown, C. T., Brooks, B., & Banfield, J. F.
>
> _The ISME Journal_, 11(12), 2864–2868. https://doi.org/10.1038/ismej.2017.126

> **Comprehensive taxonomic identification of microbial species in metagenomic data using SingleM and Sandpiper. (2025).**
>
> Woodcroft, B. J., Aroney, S. T. N., Zhao, R., Cunningham, M., Mitchell, J. A. M., Nurdiansyah, R., Blackall, L., & Tyson, G. W.
>
> _Nature Biotechnology_, 1–6. https://doi.org/10.1038/s41587-025-02738-1

> **Fast gapped-read alignment with Bowtie 2. (2012).**
>
> Langmead, B., & Salzberg, S. L.
>
> _Nature Methods_, 9(4), 357–359. https://doi.org/10.1038/nmeth.1923

> **CoverM: Read alignment statistics for metagenomics. (2025).**
>
> Aroney, S. T. N., Newell, R. J. P., Nissen, J. N., Camargo, A. P., Tyson, G. W., & Woodcroft, B. J.
>
> _Bioinformatics_, 41(4), btaf147. https://doi.org/10.1093/bioinformatics/btaf147

> **Prodigal: Prokaryotic gene recognition and translation initiation site identification. (2010).**
>
> Hyatt, D., Chen, G.-L., LoCascio, P. F., Land, M. L., Larimer, F. W., & Hauser, L. J.
>
> _BMC Bioinformatics_, 11(1), 119. https://doi.org/10.1186/1471-2105-11-119

> **MetaCerberus: Distributed highly parallelized HMM-based processing for robust functional annotation across the tree of life. (2024).**
>
> Figueroa III, J. L., Dhungel, E., Bellanger, M., Brouwer, C. R., & White III, R. A.
>
> _Bioinformatics_, 40(3), btae119. https://doi.org/10.1093/bioinformatics/btae119

> **inStrain profiles population microdiversity from metagenomic data and sensitively detects shared microbial strains. (2021).**
>
> Olm, M. R., Crits-Christoph, A., Bouma-Gregson, K., Firek, B. A., Morowitz, M. J., & Banfield, J. F.
>
> _Nature Biotechnology_, 39(6), 727–736. https://doi.org/10.1038/s41587-020-00797-0

> **Ecologically coherent population structure of uncultivated bacterioplankton. (2021).**
>
> Sjöqvist, C., Delgado, L. F., Alneberg, J., & Andersson, A. F.
>
> _The ISME Journal_, 15(10), 3034–3049. https://doi.org/10.1038/s41396-021-00985-z
