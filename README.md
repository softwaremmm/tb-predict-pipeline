[![Tests](https://github.com/oxfordmmm/tb-predict-pipeline/actions/workflows/tests.yaml/badge.svg)](https://github.com/oxfordmmm/tb-predict-pipeline/actions/workflows/tests.yaml)
# TB-Predict-Pipeline
Nextflow pipeline for producing variants, mutations and effects of a specified (minos) VCF file

## Requirements
- [Docker](https://docs.docker.com/get-docker/)
- Nextflow


## Running the NextFlow
Workflow takes parameters:
- seq_platform. `ont` or `illumina` (`illumina` by default)
- sample. path to vcf file, likely from minos
- gvcf. Path to gvcf.
- reference. reference gbk file
- catalogue. Path to mutations catalogue
- null_position. Path to list of null positions.

To save output files need to set `--publish true` which will save output files to `results`.

Example running locally:
```
nextflow run . -profile local --publish true --seq_platform illumina \
    --sample tests/test-cases/NC_045512.2-S_E484K-minos.vcf \
    --gvcf tests/test-cases/empty.gvcf \
    --reference tests/test-cases/NC_045512.2.gbk \
    --catalogue tests/test-cases/NC_045512.2-test-catalogue.csv \
    --null_positions tests/test-cases/no-null-positions.txt
```

### Batch locally
To run a batch of files locally though this pipeline requires input files to be a folder structure
- directory
    - sample1
        - minos.vcf
        - minos.gvcf
    - sample2
        - minos.vcf
        - minos.gvcf

Then use entry `batch`, and provide globs for the samples and gvcfs

```bash
nextflow run . -entry batch -profile local --publish true --seq_platform illumina \
    --samples "tests/batch_dir/*/*.vcf" \
    --gvcfs "tests/batch_dir/*/*.gvcf" \
    --reference tests/test-cases/NC_045512.2.gbk \
    --catalogue tests/test-cases/NC_045512.2-test-catalogue.csv \
    --null_positions tests/test-cases/no-null-positions.txt
```

## Tests
As this is a very simple pipeline which just calls `gnomonicus`, the majority of testing occurs within `gnomonicus`. However, the included tests cover some edge cases of files being produced, as well as ensuring files are placed in correct places.
### Run the tests
Pipelines and evaluation are run using a bash script. This should also be automatically run on each push to the `main` branch through a CI action. Note that this repo intentionally does not use `nf-test` for these tests. This is due to a few reasons, but mostly that direct file comparisons are not reliable here due to fields such as timestamp
```
tests/test.sh
```


## Run using Nextflow's kuberun
Kuberun is the only way I've found to get this running with Kubernetes. However, I have only tested this on an SP3 stack with Kubernetes setup so YMMV.
An SP3 stack has a directory `/data` which corresponds to a persistent volume claim `default-nextflow-fss-storage-data-pvc`. When Nextflow runs `kuberun`, this directory can be mounted (in this case to `/data`) on the pod - giving a working directory of `/data/<user>` (in the SP3 stack case, `/data/ubuntu`).
If sample VCFs, reference genomes and catalogues are copied/moved/linked to folders within `/data`, they can be referenced as paths within the Nextflow command:
```
nextflow kuberun https://github.com/oxfordmmm/tb-predict-pipeline -r k8sTest -v default-nextflow-fss-storage-data-pvc:/data --sample /data/<VCF path> --reference /data/<reference path> --catalogue /data/<catalogue path> --output_dir /data/<output path>
```
The result of this is utilising the specified VCF, reference and sample to populate the specified output directory within the `/data` path.

## Run with Kubernetes
This does not currently work due to unknown issues with Nextflow and Kubernetes...
```
nextflow run . -profile kubernetes --sample <absolute sample path> --reference <absolute reference path> --catalogue <absolute catalogue path> --output_dir <absolute output directory path>
```

### Kubernetes problems...
If Kubernetes can be installed and setup correctly (I had issues with installing on bare metal using 5.18.15-arch1-1), there are other issues with Nextflow integration.


Using an OCI VM and setting up with `minikube` (see below) results in error messages complaining about `.command.run not found` - presumably an issue with Nextflow not copying this to the pod.


Using a bare SP3 stack:
* Using `kubectl` at all will result in warnings about cryptographic library depreciation due to OCI using python 3.6 which reached EOL Dec 2021
* Attempting to run Nextflow using Kubernetes will result in the warnings being printed along with a malformed base-64 string (describing the pod to be created) as an error
* This can be resolved with `export PYTHONWARNINGS="ignore"` to remove this issue
* But this lands you back at the same issue of `.command.run not found`...


[This](https://github.com/seqeralabs/nf-k8s-best-practices) suggests that the `.command.run not found` is an issue with pvc latency. However, according to `kubectl describe <failed pod>`, the command being run has the path to the host's work directory (i.e it is not being copied/linked to the pod's pvc). Logging into the head pod to inspect the pvc also shows nothing being copied.


Examples of the stack traces for some issues I have found can be found in the `k8sErrorExamples` directory


#### Setting up minikube
```
curl -LO https://storage.googleapis.com/minikube/releases/latest/minikube_latest_amd64.deb
sudo dpkg -i minikube_latest_amd64.deb
sudo apt install docker.io
sudo systemctl enable --now docker
sudo usermod -aG docker $USER && newgrp docker
*** restart or log out for user changes (or use su - $USER) ***
minikube start
```

## Tags, Releases, and Committing
Use conventional commits. This is enforced with commitizen validate action and pre-commit hooks:
```bash
pre-commit install
```

This repo uses a standard gitflow approach, so changes should be first merged into develop and then released to main.
- In the develop branch semantic versioning is not used. Instead you can reference the commit hash to use it in a workflow.
- In a release branch you can create a release candidate with `cz bump a.b.c-rcX`. This also creates a tag.
- When release branch is ready for main run `cz bump a.b.c --files-only`. Manually write a human descriptive changelog. Then push these changes to main and make a release/tag there.
