import pytest
import pandas as pd
from src.summarization.gene_summarizer import rank_genes, filter_curated

def test_rank_genes():
    # Mir101c: 5 enh, 7 var, dist 3793
    # Atp5g2: 1 enh, 1 var, dist 266
    data = {
        'gene_name': ['Mir101c', 'Atp5g2'],
        'n_enhancers': [5, 1],
        'n_variants': [7, 1],
        'min_distance': [3793, 266]
    }
    df = pd.DataFrame(data)
    ranked = rank_genes(df)
    
    # Mir101c should be #1
    assert ranked.iloc[0]['gene_name'] == 'Mir101c'
    assert ranked.iloc[0]['priority_score'] > ranked.iloc[1]['priority_score']

def test_filter_curated():
    # One real gene, one Gm (predicted), one far away
    data = {
        'gene_name': ['Atp5g2', 'Gm12345', 'FarGene'],
        'distance': [1000, 1000, 100000]
    }
    df = pd.DataFrame(data)
    filtered = filter_curated(df, max_dist=50000)
    
    # Should only keep Atp5g2
    assert len(filtered) == 1
    assert filtered.iloc[0]['gene_name'] == 'Atp5g2'
