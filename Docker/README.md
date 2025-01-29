# seq_typing - Docker

**seq_typing** - Determines which reference sequence is more likely to be present in a given sample

<https://github.com/B-UMMI/seq_typing>

---

This is a Dockerfile to produce a **seq_typing** image, with all dependencies already installed.

Within this image you can find:
- ubuntu:16.04
- git v2.7.4
- Python v2.7.12 and v3.5.2
- [Blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi) v2.7.1
- [ReMatCh](https://github.com/B-UMMI/ReMatCh) v4.0.1
- [seq_typing](https://github.com/B-UMMI/seq_typing) v2.3


### Using play-with-docker
[![Try in PWD](https://cdn.rawgit.com/play-with-docker/stacks/cff22438/assets/images/button.png)](http://labs.play-with-docker.com/)

Within [play-with-docker](http://labs.play-with-docker.com/) webpage click on **create session**. Then, another page
will open with a big counter on the upper left corner. Click on **+ add new instance** and a terminal like instance should be generated on the right. On
this terminal you can load this docker image as follows:

`docker pull ummidock/seq_typing:3.0.0`

#### Build this docker on your local machine

For this, docker needs to be installed on your machine. Instructions for this can be found [here](https://docs.docker.com/engine/installation/).

##### Using DockerHub (automated build image)

`docker pull ummidock/seq_typing:3.0.0`

##### Using GitHub (build docker image)

1) `git clone https://github.com/B-UMMI/seq_typing.git`  
2) `docker build -t ummidock/seq_typing:3.0.0 ./seq_typing/Docker/`

### Run (using automated build image)
Example of _Haemophilus influenzae_ serotyping using reads and provided references sequences (adapted from [here](../README.md#reads)).
````bash
docker run --rm -u $(id -u):$(id -g) -v /local/folder/fastq_data:/data/ ummidock/seq_typing:3.0.0 \
    seq_typing.py reads --org Haemophilus influenzae \
                        --fastq /data/sample_1.fq.gz /data/sample_2.fq.gz \
                        --outdir /data/sample_out/ \
                        --threads 2
````
For more examples on how to run **seq_typing** without a Docker container see general [README](../README.md) file.


### udocker

> "A basic user tool to execute simple docker containers in user space without requiring root privileges.". From [here](https://github.com/indigo-dc/udocker).

```bash
# Get Docker image
udocker pull ummidock/seq_typing:3.0.0

# Create container (only needed to be done once)
udocker create --name=seq_typing_3-0-0 ummidock/seq_typing:3.0.0

# Run seq_typing
udocker run --user $(id -u):$(id -g) -v /local/folder/fastq_data:/data/ seq_typing_3-0-0 \
    seq_typing.py reads --org Haemophilus influenzae \
                        --fastq /data/sample_1.fq.gz /data/sample_2.fq.gz \
                        --outdir /data/sample_out/ \
                        --threads 2
```
More examples on how to use **udocker** can be found in **udocker** [GitHub page](https://github.com/indigo-dc/udocker)  
  
*__NOTE__*: if some `Error: in download: HTTP/1.1 400 Bad Request` occur while pulling the Docker image, first pull it using `docker` and then retry pulling it with `udocker` (same with Shifter).


Contact
-------
Miguel Machado <mpmachado@medicina.ulisboa.pt>  
