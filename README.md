# UDN-analysis
This repository is associated with publication "Finding Commonalities in rare diseases through the Undiagnosed Diseases Network" [ref](xxx). It allows to conduct analyses on the Undiagnosed Diseases Network (UDN) database, and cluster patients harnessing phenotypic similarity using the Louvain method. [[1]](#1).

## Prerequisites
The following libraries must be installed: 
```

argparse
sys
numpy
matplotlib
scipy
networkx
python-louvain
seaborn
collections
os
xml
operator
pandas
csv
sklearn
docx 
ast

PicSureHpdsLib
PicSureClient

```

The libraries PicSureHpdsLib and PicSureClient are from a working repository. Here are the repository links to [PicSureHpdsLib](https://github.com/hms-dbmi/pic-sure-python-adapter-hpds) and [PicSureClient](https://github.com/hms-dbmi/pic-sure-python-client)


## How to
The main body of the analysis is in the Data_analysis_UDN.py. The other scripts contain the functions necessary for computation, broken down according to the part of the analysis. The zip file xml_orphanet must be unzipped before usage, and the folder must be in the same working directory as Data_analysis_UDN.py. 

Example usage from the cmd is 
```
python Data_analysis_UDN.py --token personal_token --json_file "path/to/file" --genes_file "path/to/gene/info" 
-- variants_file "path/to/variant/info" 
```
## Authorizations
To have access to the UDN data, you will need IRB approval. When it is the case, you will need a token to log into the database, to add in the PicSureClient to connect.

## Publication
This code supports the analysis presented in "Finding Commonalities in rare diseases through the Undiagnosed Diseases Network"(publication to come).

## References
<a id="1">[1]</a> 
Blondel et al. (2008). 
Fast unfolding of communities in large networks. 
Journal of Statistical Mechanics: Theory and Experiment, Volume 2008, October 2008.

## License

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
