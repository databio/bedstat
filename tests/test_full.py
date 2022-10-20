import bedstat
import os
import pytest
import bbconf


file_path = "./tests"
# file_path = "./"
# CONFIG_FILE = f"{file_path}/config_db_github.yaml"
CONFIG_FILE = f"{file_path}/config_db_local.yaml"
SAMPLE_YAML = f"{file_path}/data/f1/AML_db358_sample_new.yaml"

# has to be change manually
OPEN_MATRIX19 = "/home/bnt4me/Virginia/bed_base_all/bedbase/bedbase_tutorial/openSignalMatrix_hg19_percentile99_01_quantNormalized_round4d.txt.gz"


@pytest.mark.skipif(True, reason="This test can be run only locally")
class TestBedstatProcessing:
    def test_processing(self):
        bedstat.run_bedstat(
            bedfile=f"{file_path}/data/f1/AML_db358.bed.gz",
            sample_yaml=SAMPLE_YAML,
            bedbase_config=CONFIG_FILE,
            bigbed=f"{file_path}/data/f1/",
            just_db_commit=False,
            genome_assembly="hg19",
            open_signal_matrix=OPEN_MATRIX19,
        )

    # We should add some additionall tests here
