import pytest
import scanpy as sc
import spatialdata as sd
import spatialdata_plot
from spatialdata import SpatialData
from os.path import join
from spatialdata_plot.pl.basic import PlotAccessor
from easy_vitessce.VitessceSpatialData import VitessceSpatialData
from easy_vitessce.configuring_plots import configure_plots

from easy_vitessce.configuring_plots import (_monkeypatch_spatialdata, _undo_monkeypatch_spatialdata)

def test_monkeypatch_sd():
    spatialdata_filepath = join("data", "merfish.spatialdata.zarr")
    sdata = sd.read_zarr(spatialdata_filepath)
    _monkeypatch_spatialdata()
    assert isinstance(sdata.pl, VitessceSpatialData)

def test_monkeypatch_sd_2():
    spatialdata_filepath = join("data", "merfish.spatialdata.zarr")
    sdata = sd.read_zarr(spatialdata_filepath)
    _undo_monkeypatch_spatialdata()
    configure_plots(enable_plots="spatial")
    assert isinstance(sdata.pl, VitessceSpatialData)    

def test_undo_monkeypatch_sd():
    spatialdata_filepath = join("data", "merfish.spatialdata.zarr")
    sdata = sd.read_zarr(spatialdata_filepath)
    _monkeypatch_spatialdata()
    print(sdata.pl)
    _undo_monkeypatch_spatialdata()
    print(sdata.pl)
    assert isinstance(sdata.pl, PlotAccessor)
    # assert 'spatialdata_plot.pl.basic.PlotAccessor' in str(sdata.pl)

def test_undo_monkeypatch_2():
    spatialdata_filepath = join("data", "merfish.spatialdata.zarr")
    sdata = sd.read_zarr(spatialdata_filepath)
    _monkeypatch_spatialdata()
    _undo_monkeypatch_spatialdata()
    assert SpatialData._orig_pl is PlotAccessor
    assert 'spatialdata_plot.pl.basic.PlotAccessor' in str(SpatialData._orig_pl)

def test_undo_monkeypatch_3():
    spatialdata_filepath = join("data", "merfish.spatialdata.zarr")
    sdata = sd.read_zarr(spatialdata_filepath)
    _monkeypatch_spatialdata()
    _undo_monkeypatch_spatialdata()
    sdata = sd.read_zarr(spatialdata_filepath) # ????
    assert isinstance(sdata.pl, PlotAccessor)

def test_embedding_config_creation():
    adata = sc.datasets.pbmc68k_reduced()
    vw = sc.pl.embedding(adata, basis="umap")
    vc = vw.config
    vc_dict = vc.to_dict(base_url='')
    assert vc_dict["datasets"][0]["files"][0]["url"].endswith(".adata.zarr")
    del vc_dict["datasets"][0]["files"][0]["url"]
    assert vc_dict == {
        'version': '1.0.15',
        'name': 'UMAP',
        'description': '',
        'datasets': [{'uid': 'A',
        'name': 'data',
        'files': [{'fileType': 'anndata.zarr',
            'options': {'obsEmbedding': [{'path': 'obsm/X_umap',
            'dims': [0, 1],
            'embeddingType': 'UMAP'}],
            'obsFeatureMatrix': {'path': 'X'}}}]}],
        'coordinationSpace': {'dataset': {'A': 'A'},
        'embeddingType': {'A': 'UMAP'},
        'featureSelection': {'A': ['']},
        'obsColorEncoding': {'A': 'geneSelection'},
        'embeddingObsRadiusMode': {'A': 'manual'},
        'embeddingObsRadius': {'A': 2.5},
        'featureValueColormap': {'A': 'viridis'}},
        'layout': [{'component': 'scatterplot',
        'coordinationScopes': {'dataset': 'A',
        'embeddingType': 'A',
        'featureSelection': 'A',
        'obsColorEncoding': 'A',
        'embeddingObsRadiusMode': 'A',
        'embeddingObsRadius': 'A',
        'featureValueColormap': 'A'},
        'x': 0,
        'y': 0,
        'w': 6,
        'h': 12},
        {'component': 'featureList',
        'coordinationScopes': {'dataset': 'A',
        'featureSelection': 'A',
        'obsColorEncoding': 'A',
        'embeddingObsRadiusMode': 'A',
        'embeddingObsRadius': 'A',
        'featureValueColormap': 'A'},
        'x': 6,
        'y': 0,
        'w': 6,
        'h': 12}],
        'initStrategy': 'auto'
    }

def test_sc_tl():
    adata = sc.datasets.pbmc68k_reduced()
    sc.pl.embedding(adata, basis="tsne")
    assert "X_tsne" in adata.obsm

@pytest.mark.skip(reason = "problem with sample_id?")
def test_spatial_config_creation():
    adata = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node", include_hires_tiff=True)
    vw = sc.pl.spatial(adata, color = "log1p_n_genes_by_counts")
    vc = vw.config
    vc_dict = vc.to_dict(base_url='')

    assert vc_dict == {
        'version': '1.0.17',
        'name': 'AnnData with image',
        'description': '',
        'datasets': [{'uid': 'A',
        'name': 'My dataset',
        'files': [{'fileType': 'anndata.zarr',
            'url': 'http://localhost:8001/A/0/28cc3238-605e-492f-8d94-eee3f88a10ca.adata.zarr',
            'options': {'obsSpots': {'path': 'obsm/spatial'},
            'obsFeatureMatrix': {'path': 'X'}}},
        {'fileType': 'image.ome-zarr',
            'url': 'http://localhost:8001/A/1/34408abc-7964-4844-bcbe-2ee24d5f00f4.ome.zarr',
            'options': {'coordinateTransformations': [{'type': 'translation',
            'translation': [0, 0, 0]},
            {'type': 'scale',
            'scale': [1, 5.878500103050106, 5.878500103050106]}]}},
        {'fileType': 'anndata.zarr',
            'url': 'http://localhost:8001/A/2/b57bcc7a-f416-4619-8296-c4f2625fe225.adata.zarr',
            'options': {'obsFeatureColumns': [{'path': 'obs/log1p_n_genes_by_counts'}]},
            'coordinationValues': {'obsType': 'cell',
            'featureType': 'qualityMetric',
            'featureValueType': 'value'}}]}],
        'coordinationSpace': {'dataset': {'A': 'A'},
        'spotLayer': {'A': '__dummy__'},
        'obsType': {'A': 'cell', 'B': 'cell'},
        'spatialSpotRadius': {'A': 45},
        'featureValueColormap': {'A': 'viridis'},
        'featureType': {'A': 'qualityMetric'},
        'featureValueType': {'A': 'value'},
        'featureSelection': {'A': ['log1p_n_genes_by_counts']},
        'obsColorEncoding': {'A': 'geneSelection'},
        'metaCoordinationScopes': {'A': {'spotLayer': ['A'],
        'obsType': 'B',
        'featureType': 'A',
        'featureValueType': 'A',
        'featureSelection': 'A',
        'obsColorEncoding': 'A'}},
        'metaCoordinationScopesBy': {'A': {'spotLayer': {'obsType': {'A': 'A'},
            'spatialSpotRadius': {'A': 'A'},
            'featureValueColormap': {'A': 'A'}}}}},
        'layout': [{'component': 'spatialBeta',
        'coordinationScopes': {'dataset': 'A',
        'metaCoordinationScopes': ['A'],
        'metaCoordinationScopesBy': ['A']},
        'x': 0,
        'y': 0,
        'w': 6,
        'h': 12},
        {'component': 'layerControllerBeta',
        'coordinationScopes': {'dataset': 'A',
        'metaCoordinationScopes': ['A'],
        'metaCoordinationScopesBy': ['A']},
        'x': 6,
        'y': 0,
        'w': 6,
        'h': 6},
        {'component': 'featureValueHistogram',
        'coordinationScopes': {'dataset': 'A'},
        'x': 6,
        'y': 6,
        'w': 6,
        'h': 6}],
        'initStrategy': 'auto'
    }