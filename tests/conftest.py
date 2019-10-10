import pytest
import os


@pytest.fixture
def data_path():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


@pytest.fixture
def csv_data_path(data_path):
    return os.path.join(data_path, "csv")


@pytest.fixture
def csv_base_dir(csv_data_path):
    return os.path.join(csv_data_path, "hg38", "test_data_csv")


@pytest.fixture
def tsv_data_path(data_path):
    return os.path.join(data_path, "tsv")


@pytest.fixture
def tsv_base_dir(tsv_data_path):
    return os.path.join(tsv_data_path, "hg38", "test_data_tsv")