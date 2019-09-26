# bedstat
EPISB pipeline for obtaining statistics about bed files

## HOW TO USE

Run "scripts/process_LOLA.py" and feed it the location of your LOLA database copy. E.g.:

```
python3 scripts/process_LOLA.py --lola_loc /ext/qumulo/LOLAweb/databases/LOLACore > lolacore.csv
```

if you pass the --genome switch, it will use that genome's folder to find all bed files. If not, it assumes "hg38".

The above command will build the csv file looper needs to run the pipeline on all the sample files from LOLA.

Then simply run the looper command to invoke the pipeline:

```
looper run project/bedstat_config.yaml
```
