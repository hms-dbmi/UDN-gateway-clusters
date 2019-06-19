# UDN-analysis
Jupyter Notebook for conducting analysis on UDN data and creation and analysis of clusters

## Missing data
The genomic data files are accessed through the JSON raw file. If you have approved access to the data, I can forward you the associated files used for genomic and variant data.
The clusters are computed with a non deterministic algorithms. To have the exact clusters analyze, if you have approved access to the data, I can forward you the associated file.

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
## Authorizations
To have access to the UDN data, you will need IRB approval. When it is the case, you will need a token to log into the database, to add in the PicSureClient to connect.

## How to
The jupyter notebook Data_analysis_UDN comprises all the code used for analyses of the UDN database.

The file clusters_un.txt is the list of clusters used for the analysis. Community detection is non-deterministic, so these clusters are susceptable to change if code is re-run. 

To use the clusters analyzed in the publication, replace ``` clusters_un ``` by that downloaded from clusters_un.txt when the community-louvain method is run. 

## Publication
This code supports the analysis presented in the following article : XXXX

## License

None for now
