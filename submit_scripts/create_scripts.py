#!/usr/bin/python

#clusterNumbers = ['600', '601', '602', '603', '604', '605', '607', '608', '609', '610']
clusterNumbers = ['602']


for cNumber in clusterNumbers:
	file = open('submit_sn_'+cNumber+'.sh',"w")
	file.write('#!/bin/bash\n')
	file.write('../sn'+ \
	' /RAID3/markus/clusterdata/galaxy_files/'+cNumber+'/new' \
	'/ /RAID3/markus/clusterdata/metal_files/'+cNumber+'/\n')
	file.close()


for cNumber in clusterNumbers:
	file = open('start_'+cNumber+'.sh',"w")
	file.write('#!/bin/bash\n')
	file.write('\n')
	file.write('#$ -q intel.q\n')
	file.write('#$ -o output_'+cNumber+'.out\n')
	file.write('#$ -cwd\n')
	file.write('#$ -j yes\n')
	file.write('./submit_sn_'+cNumber+'.sh\n')
	file.close()


#submit script
file = open('all_submit.sh','w')
file.write('#!/bin/bash\n')
for cNumber in clusterNumbers:
	file.write('qsub start_'+cNumber+'.sh\n')
	file.write('sleep 0.5\n')
file.close()
