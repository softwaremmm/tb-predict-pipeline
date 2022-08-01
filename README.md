# TB-Predict-Pipeline
Nextflow pipeline for producing variants, mutations and effects of a specified (minos) VCF file

## Run locally with Docker
```
nextflow run . -profile docker --sample <absolute sample path> --reference <absolute reference path> --catalogue <absolute catalogue path> --output_dir <absolute output directory path>
```

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