# MZML Converter
https://mzcloudconv.uc.r.appspot.com/

Lightweight web app built in Flask and deployed on Google Cloud. Uses the [MassBank NIST dataset](https://github.com/MassBank/MassBank-data) for compound reference.

## Structure
`app.py`- Flask backend and GCS storage. All uploaded files are stored in Google Cloud.

`data_parser.py`- Module for loading the MNIST dataset and parsing the uploaded mzML into a new XML file.

`predictor.py`- Module that identifies the peaks in the MS scans and finds the closest matching compound. Uses a custom cosine similarity algorithm for predictions. Also extracts instrument information.

## mzML files used to test + Disclaimer
Smaller mzML files were used to test and develop this app.

[01R1.mzML](https://massive.ucsd.edu/ProteoSAFe/dataset_files.jsp?task=ed793078c76a495cb0c9c2c7b7bba84b#%7B%22table_sort_history%22%3A%22main.collection_asc%22%7D),
[D1a.mzML](https://massive.ucsd.edu/ProteoSAFe/dataset_files.jsp?task=a6c218ec9c784f93b08edd1fd32ee78d#%7B%22table_sort_history%22%3A%22main.collection_asc%22%7D)
from the UCSD database are two sample files that were tested and worked very well. Very large mzML files or ones that to not strictly conform to the XML standard may not be properly parsed by the app.

Built by Benson Zhang.