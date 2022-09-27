from scripts.process_LOLA import process_index_file


def test_process_index_file(capsys, tsv_base_dir, csv_base_dir):
    """Check if process_index_file produces the same results for TSV and CSV files"""
    process_index_file(csv_base_dir, "hg38")
    out1, _ = capsys.readouterr()
    process_index_file(tsv_base_dir, "hg38")
    out2, _ = capsys.readouterr()
    # number of characters in each captured output should be exactly the same
    assert len(out1) == len(out2)
