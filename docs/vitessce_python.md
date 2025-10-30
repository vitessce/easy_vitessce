# Advanced

## Access the Vitessce configuration

For advanced usage, you can obtain the [VitessceConfig](https://python-docs.vitessce.io/api_config.html#vitessce.config.VitessceConfig) instance.

```py
sc.pl.dotplot(adata, ...)
```

First, modify the plotting function call so that the [VitessceWidget](https://python-docs.vitessce.io/api_config.html#vitessce.widget.VitessceWidget) instance can be accessed via a variable:

```py
vw = sc.pl.dotplot(adata, ...)
vw
```

The `vw` variable contains the `VitessceWidget` instance, which allows inspecting the `VitessceConfig` instance.

In a subsequent notebook cell:

```py
# Access the configuration powering this widget instance.
vw.config
# or
vw.config.to_dict(base_url="")
# or, for "live" dynamically-updated config
vw._config
```

### Access values from the coordination space

For background, see the documentation on [coordinated multiple views in Vitessce](https://vitessce.io/docs/coordination/).

For instance, to access the list of currently-selected genes (i.e., features):

```py
coordination_type = "featureSelection"
coordination_object = vw._config["coordinationSpace"][coordination_type]
coordination_object
```

Or, to access the set of currently-selected cell types (i.e., observation sets):

```py
coordination_type = "obsSetSelection"
coordination_object = vw._config["coordinationSpace"][coordination_type]
coordination_object
```

See the full list of [coordination types](https://vitessce.io/docs/coordination-types/) for more options.

## When to use `vitessce` over `easy_vitessce`

An advantage of using the plain `vitessce` Python [package](https://python-docs.vitessce.io/) is that you do not need to have your AnnData/SpatialData object present/accessible.
In fact, you can configure vitessce for hypothetical objects that do not exist or will never exist.
This is useful if you want to construct configurations for objects that will be materialized later, e.g., in a data portal context.
However, it also has the downside that you can easily construct an erroneous config (e.g., you specify a column name that does not exist in the targeted AnnData object) as the error will not arise until runtime/rendering-time.

On the other hand, `easy_vitessce` uses the object contents for certain logic.
For example, when you specify `color="GeneA"` vs. `color="cell_type"` vs. `color="red"` in plotting params, we perform checks like:

```py
if color in adata.var.index:
  # color was a gene name
  pass
elif color in adata.obs.columns:
  # color was an obs column name
  pass
else:
  # color may have been a color string or hex code string
  pass
```