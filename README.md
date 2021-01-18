<div align="center">

[![Logo](logo/multilayer2.png)](https://www.sysfate.org/)
</div>

<a href="https://www.biorxiv.org/content/10.1101/2020.11.09.374660v1"><img src="logo/doi_multilayer.png" width="25%" height="25%"></a>

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

## Description

Inspired by contextual pixel classification strategies applied to image analysis, we have developed MULTILAYER, 
allowing to stratify spatially-resolved transcriptome maps into functionally-relevant molecular substructures. 
For this, MULTILAYER applies agglomerative clustering within contiguous locally-defined transcriptomes (herein 
defined as gene expression elements or gexels), combined with community detection methods for graph partitioning.

## Run

Launch : `python3 Multilayer.py`
 
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

###### - Louvain

```bash
pip install python-louvain
```

###### - Python Imaging Library / PIL

```bash
pip install pillow
```

# MULTILAYER COMPRESSOR

## Description

We developed module for compress data. Multilayer compressor is able to merge several gexels in one big gexel.

## Run

Launch : `python3 Multilayer_Compressor.py`

## Dependencies

###### - Numpy

```bash
pip install numpy
```

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
