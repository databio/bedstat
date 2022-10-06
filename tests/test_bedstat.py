import bedstat
import os
import pytest
import bbconf


file_path = "./tests"
# file_path = "./"
CONFIG_FILE = f"{file_path}/config_db_github.yaml"
# CONFIG_FILE = f"{file_path}/config_db_local.yaml"


class TestBedstat:
    def test_uploading(self):
        bedstat.run_bedstat(
            bedfile=f"{file_path}/data/f1/AML_db358.bed.gz",
            sample_yaml=f"{file_path}/data/f1/AML_db358_sample.yaml",
            bedbase_config=CONFIG_FILE,
            bigbed=f"{file_path}/data/f1/",
            just_db_commit=True,
            genome_assembly="hg19",
        )

    def test_db_content(self):
        bbc = bbconf.BedBaseConf(config_path=CONFIG_FILE, database_only=True)
        assert bbc.bed.record_count == 1

    def test_db_record(self):
        bbc = bbconf.BedBaseConf(config_path=CONFIG_FILE, database_only=True)
        print(bbc.bed.select(columns=["other"]))
        geonome_check = bbc.bed.select(columns=["other"])[0]["genome"]
        assert geonome_check == "hg19"
