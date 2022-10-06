import bedstat
import os
import pytest


file_path = "./tests"
# file_path = "./"
CONFIG_FILE = f"{file_path}/config_db_github.yaml"


class TestBedstat:
    def test_uploading(self):
        bedstat.run_bedstat(
            bedfile=f"{file_path}/data/f1/AML_db358.bed.gz",
            genome_assembly="hg19",
            sample_yaml=f"{file_path}/data/f1/AML_db358_sample.yaml",
            bedbase_config=CONFIG_FILE,
            bigbed=f"{file_path}/data/f1/",
            just_db_commit=True,
        )
