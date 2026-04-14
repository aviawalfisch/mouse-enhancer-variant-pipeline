import pytest
import pandas as pd
from src.annotation.gene_mapper import map_and_score, attach_to_variants

def test_map_and_score(tmp_path):
    # Mock bedtools closest output
    # enh_chr, enh_start, enh_end, enhancer_id, gene_chr, gene_start, gene_end, gene_name, distance
    nearest_data = "chr1\t1000\t2000\tenh1\tchr1\t500\t501\tGeneA\t500\n"
    nearest_file = tmp_path / "nearest.tsv"
    nearest_file.write_text(nearest_data)
    
    # Test distance-only scoring
    df = map_and_score(str(nearest_file))
    
    assert len(df) == 1
    assert df.iloc[0]['gene_name'] == "GeneA"
    assert df.iloc[0]['score'] == 1 / 500

def test_attach_to_variants(tmp_path):
    # Scored enhancers
    scored_data = {
        'enh_chr': ['chr1'],
        'enh_start': [1000],
        'enh_end': [2000],
        'gene_name': ['GeneA'],
        'score': [0.1]
    }
    scored_df = pd.DataFrame(scored_data)
    
    # Variants file
    variants_data = "variant1\tchr1\t1500\tA\tG\t.\t.\t.\tchr1\t1000\t2000\tannot1\n"
    # Column order from variant_parser: variant_id, chr, pos, ref, alt, qual, filter, info, enh_chr, enh_start, enh_end, annotation
    cols = ["variant_id", "chr", "pos", "ref", "alt", "qual", "filter", "info", "enh_chr", "enh_start", "enh_end", "annotation"]
    variants_file = tmp_path / "variants.tsv"
    variants_file.write_text("\t".join(cols) + "\n" + variants_data)
    
    final_df = attach_to_variants(str(variants_file), scored_df)
    
    assert len(final_df) == 1
    assert final_df.iloc[0]['gene_name'] == "GeneA"
    assert final_df.iloc[0]['variant_id'] == "variant1"
