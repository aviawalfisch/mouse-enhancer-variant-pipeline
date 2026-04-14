import pytest
import pandas as pd
import os
from src.annotation.tss_generator import generate_tss_bed

def test_generate_tss_bed(tmp_path):
    # Create a dummy GTF file
    gtf_content = [
        'chr1\tHAVANA\tgene\t100\t200\t.\t+\t.\tgene_id "G1"; gene_name "GeneA";',
        'chr1\tHAVANA\tgene\t300\t400\t.\t-\t.\tgene_id "G2"; gene_name "GeneB";',
    ]
    gtf_file = tmp_path / "test.gtf"
    gtf_file.write_text("\n".join(gtf_content))
    
    df = generate_tss_bed(str(gtf_file))
    
    # Check results
    assert len(df) == 2
    # GeneA: + strand, start 100 -> TSS 99
    assert df.iloc[0]['start'] == 99
    assert df.iloc[0]['gene_name'] == "GeneA"
    
    # GeneB: - strand, end 400 -> TSS 399
    assert df.iloc[1]['start'] == 399
    assert df.iloc[1]['gene_name'] == "GeneB"

