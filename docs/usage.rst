Usage
=====

.. _installation:

Installation
------------

This document covers the installation requirements and dependency management for the PopMAG pipeline. It explains how to set up the required software environment, configure container systems, and manage external databases needed for population genomics analysis based on metagenome-assembled genomes.

For basic usage instructions after installation, see Basic Usage. For detailed configuration of pipeline parameters, see Global Parameters.

Prerequisites
-------------

* Nextflow ≥24.04.2 nextflow.config:161
* Java ≥11 (required by Nextflow)
* Container Engine: Docker, Singularity, Podman, or Apptainer
* Conda/Mamba (optional alternative to containers)

Installation Steps
------------------

#. Install Nextflow

.. code-block:: console

   curl -s https://get.nextflow.io | bash  
   chmod +x nextflow  
   sudo mv nextflow /usr/local/bin/

Verify installation: 

.. code-block:: console

    nextflow info

#. Choose Execution Profile
PopMAG provides multiple execution profiles for different environments:

* Docker (recommended): -profile docker 
* Singularity (HPC systems): -profile singularity
* Apptainer: -profile apptainer
* Mamba: -profile mamba
* Conda: -profile conda 
