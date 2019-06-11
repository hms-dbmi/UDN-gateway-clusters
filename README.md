# UDN-analysis
Jupyter Notebook for conducting analysis on UDN data and creation and analysis of clusters

## Prerequisites
The following libraries must be installed: 
```bash
numpy
scipy
matplotlib
seaborn
community
networkx
collection
os
xml.etree.ElementTree
PicSureHpdsLib
PicSureClient
```

Here are the repo links to [PicSureHpdsLib](https://github.com/hms-dbmi/pic-sure-python-adapter-hpds) and [PicSureClient](https://github.com/hms-dbmi/pic-sure-python-client)

## How to
The jupyter notebook Data_analysis_UDN comprises all the code used for analyses of the UDN database.

The file clusters_un.txt is the list of clusters used for the analysis. Community detection is non-deterministic, so these clusters are susceptable to change if code is re-run. 

To use the clusters analyzed in the publication, replace ``` clusters_un ``` by that downloaded from clusters_un.txt when the community-louvain method is run. 

## Publication
This code supports the analysis presented in the following article : XXXX

## License

None for now
