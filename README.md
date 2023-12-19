Archived: This pipeline is now maintained as part of https://github.com/databio/bedboss

# bedstat
pipeline for obtaining statistics about bed files

## Installation instructions

1. Clone this repository
2. Install the python packages listed in the `requirements.txt` file

```
pip install -r requirements.txt --user
```
3. Install additional dependacies listed below

## How to use

### 0. (optional) Convert a LOLA database into a PEP

*This step is required only if you start from a "LOLA Region Database"*

Run `scripts/process_LOLA.py` and feed it the location of your LOLA database copy. E.g.:

```
python3 scripts/process_LOLA.py --lola_loc /ext/qumulo/LOLAweb/databases/LOLACore > lolacore.csv
```

if you pass the `--genome` switch, it will use that genome's folder to find all bed files. If not, it assumes "hg38".

The above command will build the csv file looper needs to run the pipeline on all the sample files from LOLA.

### 1. Validate your PEP with [`eido`](https://github.com/pepkit/eido)

The input PEP can be validated against the [JSON schema in this repository](pep_schema.yaml). This ensures the PEP consists of all required attributes to run `bedstat` pipeline.

```
eido validate <path/to/pep> -s https://schema.databio.org/pipelines/bedstat.yaml
```

### 2. Create a persistent volume to house elasticsearch data

```
docker volume create es-data
```

### 3. Run the docker container for elasticsearch

```
docker run -p 9200:9200 -p 9300:9300 -v es-data:/usr/share/elasticsearch/data -e "xpack.ml.enabled=false" \
  -e "discovery.type=single-node" elasticsearch:7.5.1
```

### 4. Run the bedstat pipeline on the PEP

Then simply run the looper command to run the pipeline for each bed file. It will create a set of plots and statistics per bed file and insert the metadata into Elastic:

```
looper run project/bedstat_config.yaml
```

The data loaded into elasticsearch should persist between elasticsearch invocations, on the es-data docker volume created above in step 2.

### 5. (optional) Run Kibana

Kibana can be used in order to see ElasticSearch data in a "GUI" kind of a way.

Pull a matching Kibana docker image. Make sure the Elasticsearch and Kibana container tags match:
```
docker pull docker.elastic.co/kibana/kibana:7.5.1
```

Get the ID of the docker container (started above) running ElasticSearch via 
```
docker ps | grep elasticsearch
```

Run Kibana to link to that container:
```
docker run --link <ID OF ELASTIC CONTAINER HERE>:elasticsearch -p 5601:5601  docker.elastic.co/kibana/kibana:7.5.1
```

Point your local web browser to http://localhost:5601

---

## Additional dependencies

[regionstat.R](tools/regionstat.R) script is used to calculate the bed file statistics, so the pipeline also depends on several R packages:

* `BiocManager`
* `optparse`
* `devtools`
* `GenomicRanges`
* `GenomicDistributions`
* `BSgenome.<organim>.UCSC.<genome>` *depending on the genome used* 
* `LOLA`

you can use [installRdeps.R](scripts/installRdeps.R) helper script to easily install the required packages:

```
Rscript scripts/installRdeps.R
``` 

