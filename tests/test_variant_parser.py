import pytest
import pandas as pd
from src.parsing.variant_parser import parse_enhancer_hits

def test_parse_enhancer_hits(tmp_path):
    # Mock input: 9 VCF columns + 6 sample columns + 4 annotation columns = 19 columns
    vcf_base = ["chr1", "1000", "rs1", "A", "G", ".", "PASS", ".", "GT"]
    samples = ["0/1"] * 6
    annot = ["chr1", "900", "1100", "enh1"]
    line = "\t".join(vcf_base + samples + annot)
    
    input_file = tmp_path / "hits.txt"
    input_file.write_text(line)
    
    v_df, e_df = parse_enhancer_hits(str(input_file), n_samples=6)
    
    # Check variants
    assert len(v_df) == 1
    assert v_df.iloc[0]['variant_id'] == "1:1000:A>G"
    assert v_df.iloc[0]['enh_start'] == 900
    
    # Check enhancers
    assert len(e_df) == 1
    assert e_df.iloc[0]['enhancer_id'] == "chr1:900-1100"
