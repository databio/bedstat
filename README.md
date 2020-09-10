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

### 2. Run PostgreSQL

For example, to run an instance in a container and make the data persist, execute:

```
docker volume create postgres-data
docker run -d --name bedbase-postgres -p 5432:5432 -e POSTGRES_PASSWORD=bedbasepassword -e POSTGRES_USER=postgres -e POSTGRES_DB=postgres -v postgres-data:/var/lib/postgresql/data postgres
```
Provided environment variables need to match the settings in bedbase configuration file

### 3. Run the bedstat pipeline on the PEP

Then simply run the `looper run` command to run the pipeline for each bed file. It will create a set of plots and statistics per bed file and insert the metadata into PostgreSQL:

```
looper run project/bedstat_config.yaml
```

The data loaded into PostgreSQL should persist between PostgreSQL invocations, on the `postgres-data` docker volume created above in step 2.

## Additional dependencies

[regionstat.R](tools/regionstat.R) script is used to calculate the bed file statistics, so the pipeline also depends on several R packages:

* `R.utils`
* `BiocManager`
* `optparse`
* `devtools`
* `GenomicRanges`
* `GenomicFeatures`
* `ensembldb`
* `GenomicDistributions`
* `BSgenome.<organim>.UCSC.<genome>` *depending on the genome used* 
* `LOLA`

you can use [installRdeps.R](scripts/installRdeps.R) helper script to easily install the required packages:

```
Rscript scripts/installRdeps.R
``` 

