#!/usr/bin/python

import os

#clusterNumbers = ['600', '601', '602', '603', '604', '605', '607', '608', '609', '610']
clusterNumbers = ['602']

for cNumber in clusterNumbers:
	path='/RAID3/markus/clusterdata/metal_files/'+cNumber+'/'
	os.mkdir( path);
	print 'created '+path+'\n'

