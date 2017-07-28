# parallel_meta

PREREQUISITE
Python 2.7
pip



INSTALLATION
$ unzip parallel_meta.zip
$ cd parallel_meta
$ pip install -r requirements.txt


SAMPLE FORMAT REQUIREMENT
under your sample folder, there should be 
	1) “meta.txt” containing meta information of your sample
	2) a subfolder named “samples” under which are 
		*) each sample should have a separate sub-folder named after the sample sequence ID
			*) within the subfolder, a “classification.txt” containing raw data of the sample
	3) a “database.txt”  that you should copy from ../ to the current sample folder.


USAGE
$ python ./parallel_meta.py /Directory/Of/Your/Sample
