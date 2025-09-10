Usage
=====

.. _installation:

Installation
------------

This document covers the installation requirements and dependency management for the PopMAG pipeline. It explains how to set up the required software environment, configure container systems, and manage external databases needed for population genomics analysis based on metagenome-assembled genomes.

For basic usage instructions after installation, see Basic Usage. For detailed configuration of pipeline parameters, see Global Parameters.

Prerequisites
~~~~~~~~~~~~~

* Nextflow ≥ 24.04.2
* Java ≥ 17
* Bash ≥ 3.2
* Container Engine: Docker, Singularity, Podman, or Apptainer
* Conda/Mamba (optional alternative to containers)

Installation Steps
~~~~~~~~~~~~~~~~~~ 

Install Nextflow
^^^^^^^^^^^^^^^^^^^

Nextflow official installation instructions can be found in https://www.nextflow.io/docs/latest/install.html

**Self-Install**

.. code-block:: bash

   curl -s https://get.nextflow.io | bash  
   chmod +x nextflow  
   sudo mv nextflow /usr/local/bin/

**Mamba/Conda**

.. code-block:: bash

   mamba create --name nf-env bioconda::nextflow
   source activate nf-env   
 
**Verify installation** 

.. code-block:: bash

    nextflow info

**2. Choose Execution Profile**

PopMAG supports multiple containerization and environment management systems for reproducible execution. The pipeline provides several execution profiles that configure different container engines and dependency management approaches:

* Docker (recommended): -profile docker 
* Singularity (HPC systems): -profile singularity
* Apptainer: -profile apptainer
* Mamba: -profile mamba
* Conda: -profile conda 
