import pytest
import pyBigWig
import pandas as pd
import os
from src.utils.bigwig_utils import calculate_bigwig_mean

def test_calculate_bigwig_mean(tmp_path):
    # Create a small BigWig file
    bw_path = str(tmp_path / "test.bw")
    bw = pyBigWig.open(bw_path, "w")
    bw.addHeader([("chr1", 1000)])
    # Add some signal: 0-100 has 10.0, 100-200 has 20.0
    bw.addEntries(["chr1"], [0], ends=[100], values=[10.0])
    bw.addEntries(["chr1"], [100], ends=[200], values=[20.0])
    bw.close()
    
    # Create a BED file
    bed_path = str(tmp_path / "test.bed")
    with open(bed_path, "w") as f:
        f.write("chr1\t50\t150\tenh1\n") # Mean should be (50*10 + 50*20)/100 = 15.0
        
    df = calculate_bigwig_mean(bw_path, bed_path)
    
    assert len(df) == 1
    assert df.iloc[0]['enhancer_id'] == "enh1"
    assert df.iloc[0]['mean_signal'] == 15.0
