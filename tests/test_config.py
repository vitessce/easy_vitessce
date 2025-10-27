import pytest
import scanpy as sc
import spatialdata as sd
import spatialdata_plot
from spatialdata import SpatialData
from os.path import join, abspath, dirname
from spatialdata_plot.pl.basic import PlotAccessor

from vitessce import VitessceConfig
from easy_vitessce.spatialdata_plot import VitesscePlotAccessor
from easy_vitessce.configure_plots import (
    configure_plots,
    _monkeypatch_spatialdata,
    _undo_monkeypatch_spatialdata
)

# We set the matplotlib backend to something non-interactive during tests.
import matplotlib
matplotlib.use('Agg')

TEST_DIR = dirname(abspath(__file__))

def test_monkeypatch_sd():
    spatialdata_filepath = join(TEST_DIR, "data", "merfish.spatialdata.zarr")
    sdata = sd.read_zarr(spatialdata_filepath)
    _monkeypatch_spatialdata()
    assert isinstance(sdata.pl, VitesscePlotAccessor)
    _undo_monkeypatch_spatialdata()


def test_monkeypatch_sd_2():
    spatialdata_filepath = join(TEST_DIR, "data", "merfish.spatialdata.zarr")
    sdata = sd.read_zarr(spatialdata_filepath)
    _undo_monkeypatch_spatialdata()
    configure_plots(enable_plots=["spatialdata-plot"])
    assert isinstance(sdata.pl, VitesscePlotAccessor)
    assert sdata.pl is not PlotAccessor


def test_undo_monkeypatch_sd():
    spatialdata_filepath = join(TEST_DIR, "data", "merfish.spatialdata.zarr")
    sdata = sd.read_zarr(spatialdata_filepath)
    _monkeypatch_spatialdata()
    assert isinstance(sdata.pl, VitesscePlotAccessor)
    assert sdata.pl is not PlotAccessor


    # Check that creating a plot after monkeypatching returns something interactive.
    vw = sdata.pl.render_images(element="rasterized").pl.show()
    assert isinstance(vw.config, VitessceConfig)

    _undo_monkeypatch_spatialdata()

    # Undoing does not seem to affect the existing instances of the class,
    # which we verify below.
    assert isinstance(sdata.pl, VitesscePlotAccessor)
    assert sdata.pl is not PlotAccessor
    
    # However, our ._is_enabled workaround should work as expected:
    # Check that creating a plot after undoing the monkeypatch returns something static, rather than interactive.
    ax = sdata.pl.render_images(element="rasterized").pl.show(return_ax=True)
    assert isinstance(ax, matplotlib.axes.Axes)

    # Undoing should affect a new instance of the class, however.
    new_sdata = sd.read_zarr(spatialdata_filepath)
    assert not isinstance(new_sdata.pl, VitesscePlotAccessor)
    assert isinstance(new_sdata.pl, PlotAccessor)
    
    # Check that creating a plot after undoing the monkeypatch still works.
    new_ax = new_sdata.pl.render_images(element="rasterized").pl.show(return_ax=True)
    assert isinstance(new_ax, matplotlib.axes.Axes)

def test_undo_monkeypatch_2():
    # Test the class itself, rather than an sdata instance.
    _monkeypatch_spatialdata()
    assert SpatialData.pl is VitesscePlotAccessor
    assert SpatialData._orig_pl is PlotAccessor
    _undo_monkeypatch_spatialdata()
    assert SpatialData.pl is PlotAccessor
    assert not hasattr(SpatialData, "_orig_pl")


def test_embedding_config_creation():
    adata = sc.datasets.pbmc68k_reduced()
    vw = sc.pl.embedding(adata, basis="umap")
    vc = vw.config
    vc_dict = vc.to_dict(base_url='')
    assert vc_dict["datasets"][0]["files"][0]["url"].endswith(".adata.zarr")
    del vc_dict["datasets"][0]["files"][0]["url"]
    assert vc_dict == {'version': '1.0.18',
    'name': 'sc.pl.embedding for UMAP',
    'description': '',
    'datasets': [{'uid': 'A',
    'name': 'embedding data',
    'files': [{'fileType': 'anndata.zarr',
        'options': {'obsEmbedding': [{'path': 'obsm/X_umap',
            'dims': [0, 1],
            'embeddingType': 'UMAP'}]}}]}],
    'coordinationSpace': {'dataset': {'A': 'A'},
    'embeddingType': {'A': 'UMAP'},
    'embeddingObsRadiusMode': {'A': 'manual'},
    'embeddingObsRadius': {'A': 2.5}},
    'layout': [{'component': 'scatterplot',
    'coordinationScopes': {'dataset': 'A',
        'embeddingType': 'A',
        'embeddingObsRadiusMode': 'A',
        'embeddingObsRadius': 'A'},
    'props': {},
    'x': 0,
    'y': 0,
    'w': 12,
    'h': 12}],
    'initStrategy': 'auto'}

def test_sc_tl():
    adata = sc.datasets.pbmc68k_reduced()
    sc.tl.tsne(adata, random_state=1)
    sc.pl.embedding(adata, basis="tsne")
    assert "X_tsne" in adata.obsm

@pytest.mark.skip(reason = "problem with V1_Human_Lymph_Node dataset")
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

def test_spatialdata_config_creation():
    spatialdata_filepath = join(TEST_DIR, "data", "merfish.spatialdata.zarr")
    sdata = sd.read_zarr(spatialdata_filepath)
    _monkeypatch_spatialdata()

    # Check that .pl.render_something can be chained together.
    vw = sdata.pl.render_images(element="rasterized").pl.render_shapes(element="cells").pl.show()
    vc = vw.config
    vc_dict = vc.to_dict(base_url='')

    assert vc_dict["datasets"][0]["files"][0]["url"].endswith(".sdata.zarr")
    assert len(vc_dict["datasets"][0]["files"]) == 2
    
    del vc_dict["datasets"][0]["files"][0]["url"]
    assert vc_dict["datasets"][0]["files"][0] == {
        'fileType': 'spatialdata.zarr',
        'options': {
            'image': {'path': 'images/rasterized'}
        },
        'coordinationValues': {'fileUid': 'image_rasterized'}
    }

    del vc_dict["datasets"][0]["files"][1]["url"]
    assert vc_dict["datasets"][0]["files"][1] == {
        'fileType': 'spatialdata.zarr',
        'options': {
            'obsFeatureMatrix': {'path': 'tables/table/X'},
            'obsSpots': {'path': 'shapes/cells', 'tablePath': 'tables/table'}
        },
        'coordinationValues': {'obsType': 'spot', 'featureType': 'gene'}
    }
