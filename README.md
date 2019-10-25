# bedstat
EPISB pipeline for obtaining statistics about bed files

## How to use

### 1. Convert a LOLA database into a PEP

Run "scripts/process_LOLA.py" and feed it the location of your LOLA database copy. E.g.:

```
python3 scripts/process_LOLA.py --lola_loc /ext/qumulo/LOLAweb/databases/LOLACore > lolacore.csv
```

if you pass the --genome switch, it will use that genome's folder to find all bed files. If not, it assumes "hg38".

The above command will build the csv file looper needs to run the pipeline on all the sample files from LOLA.

### 2. Create a persistent volume to house elasticsearch data

```
docker volume create es-data
```

### 3. Run the docker container for elasticsearch

```
docker run -p 9200:9200 -p 9300:9300 -v es-data:/usr/share/elasticsearch/data elasticsearch:6.5.4
```

### 4. Run the bedstat pipeline on the PEP

Then simply run the looper command to invoke the pipeline:

```
looper run project/bedstat_config.yaml
```

The data loaded into elasticsearch should persist between elasticsearch invocations, on the es-data docker volume created above in step 2.

### Optional step = run Kibana

Kibana can be used in order to see ElasticSearch data in a "GUI" kind of a way.

Pull the Kibana docker image:
```
docker pull docker.elastic.co/kibana/kibana:6.5.4
```

Get the ID of the docker container (started above) running ElasticSearch via 
```
docker ps | grep elasticsearch
```

Run Kibana to link to that container:
```
docker run --link <ID OF ELASTIC CONTAINER HERE>:elasticsearch -p 5601:5601  docker.elastic.co/kibana/kibana:6.5.4
```

Point your local web browser to http://localhost:5601

### R Dependencies ###

Following R packages are necessary to run the code that processes BED files:

* BiocManager
* optparse
* devtools
* GenomicRanges (via BiocManager::install)
* GenomicDistributions (via devtools::install_github("databio/GenomicDistributions")
* BSgenome (via BiocManager::install)
* LOLA (via BiocManager::install)
