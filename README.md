<div align="center">

[![Logo](logo/multilayer2.png)](https://www.sysfate.org/)
</div>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6365536.svg)](https://www.cell.com/cell-systems/fulltext/S2405-4712(21)00151-4?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2405471221001514%3Fshowall%3Dtrue) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6365580.svg)](https://star-protocols.cell.com/protocols/998)

## Informations / Authors

|         |                                                                                               |
| ------- | --------------------------------------------------------------------------------------------- |
| Author  | MOEHLIN Julien ([Github](https://github.com/JulienMoehlin), [Gitlab](https://gitlab.com/julienmoehlin)) |
| Author  | MOLLET Bastien                                                                                |
| Author  | COLOMBO Bruno Maria                                                                           |
| Author  | MENDOZA PARRA Marco ([Github](https://github.com/SysFate))                                    |
| Team    | [SysFate](https://www.sysfate.org/)                                                           |
| Email   | <mmendoza@genoscope.cns.fr>                                                                   |


# MULTILAYER TOOL

## Changelog 2024

- The tool now detect automaticaly if the data are interstitial
- Fixed pattern detection in case of interstitial data
- Added a script to increase the number of pixels
- Improved several graphical elements

## Description

Inspired by contextual pixel classification strategies applied to image analysis, we have developed MULTILAYER, 
allowing to stratify spatially-resolved transcriptome maps into functionally-relevant molecular substructures. 
For this, MULTILAYER applies agglomerative clustering within contiguous locally-defined transcriptomes (herein 
defined as gene expression elements or gexels), combined with community detection methods for graph partitioning.

## Run

Launch : 
```python
python3 Multilayer.py
```
 
## Dependencies

###### - Numpy

```bash
pip install numpy
```

###### - Matplotlib

```bash
pip install matplotlib
```

###### - Pandas

```bash
pip install pandas
```

###### - Scipy

```bash
pip install scipy
```

###### - Scikit-learn

```bash
pip install scikit-learn
```

###### - Seaborn

```bash
pip install seaborn
```

###### - Networkx

```bash
pip install networkx
```

###### - Louvain (from Community)

```bash
pip install python-louvain
```

###### - Python Imaging Library / PIL

```bash
pip install pillow
```

#### Only needed for UMAP

###### - Louvain (from Louvain)

```bash
pip install louvain
```

###### - Scanpy

```bash
pip install scanpy
```
if you use conda :

```bash
conda install -c conda-forge scanpy python-igraph leidenalg
```

# MULTILAYER COMPRESSOR

## Description

We developed module for compress data. Multilayer compressor is able to merge several gexels in one big gexel.

## Run

Launch : 
```python
python3 Multilayer_Compressor.py -i input.tsv -o output.tsv -cx 100 -cy 100
```

## Dependencies

###### - Numpy

```bash
pip install numpy
```

###### - Pandas

```bash
pip install pandas
```

# MULTILAYER EXPANDER

## Description

We developed module for extend data. Multilayer expander is able to divide one gexel into several gexels. n is the number of expansion.

## Run

Launch : 
```python
python3 multilayer-expander.py -n 1
```

## Dependencies

###### - Pandas

```bash
pip install pandas
```


# VISIUM CONVERTER

## Description

You need cellranger of [10xgenomics](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/6.0), as well as the datasets provided by the Visium platform; namely the matrix in h5 format (Feature / cell matrix HDF5), the feature information files (Feature / cell matrix), the spatial information data (Spatial imaging data), as well as our specialized script “visium converter.py” (dependencies: pandas package).

## Run

First, use this command on “cellranger” to convert the h5 matrix to csv format: 
`./bin/cellranger mat2csv Feature/cell matrix_HDF5.h5 out_file_matrix.csv`

Then, use our python script “visiumConverter.py” as following:

```python
python3 visiumConverter.py -m out_file_matrix.csv -p spatial/tissue_positions_list.csv -g raw_feature_bc_matrix/features.tsv.gz -o matrix_multilayer.tsv –compressor
```

## Dependencies

###### - Pandas

```bash
pip install pandas
```

# ENRICHR CONVERTER

## Description

A converter for [Enrichr libraries](https://maayanlab.cloud/Enrichr/#libraries). Once the library converted, you have to place it in the directory called 'GO DB'.

## Run

Launch : 

```python
python3 enrichr_converter.py -i input.tsv
```

## Dependencies

###### - Pandas

```bash
pip install pandas
```

# DATA

All data available : (Article - data)

<ul>
<li> 

[Berglund & al](https://www.nature.com/articles/s41467-018-04724-5) - [Prostate cancer data](Data/Prostate_cancer) 
</li>

<li>

[Asp & al](https://www.sciencedirect.com/science/article/abs/pii/S0092867419312826?via%3Dihub) - [Heart development data](Data/Development_heart) 
</li>

<li>

[Moncada & al](https://www.nature.com/articles/s41587-019-0392-8) - [Pancreatic adenocarcinoma](Data/Pancreatic_adenocarcinoma) 
</li>

<li>

[Rodriques & al](https://science.sciencemag.org/content/363/6434/1463) - [High resolution brain](Data/High_resolution_brain) 
</li>

<li>

[Ståhl & al](https://science.sciencemag.org/content/353/6294/78) - [Breast cancer](Data/Breast_cancer) 
</li>

<li>

[Ståhl & al](https://science.sciencemag.org/content/353/6294/78) - [Mouse olfactory bulb](Data/Mouse_Olfactory_Bulb) 
</li>

<li>

[Liu & al](https://www.sciencedirect.com/science/article/abs/pii/S0092867420313908) - [Whole mouse embryo](Data/Whole_mouse_embryo) 
</li>
</ul>

# TUTORIAL

Multilayer [tutorial](/MULTILAYER-Tutorial.pdf)

<a href="https://www.youtube.com/channel/UCoGGql3etFGqHifFf4uDSYQ"><img src="logo/youtube_logo.png" width="30%" height="30%"></a>
