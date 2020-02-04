from multiprocessing import Process
from random import randint
from subprocess import check_output
import sys

def sample_function(rscript_filename, directory1, directory2, name, pop, size):
	check_output(["Rscript",rscript_filename, directory1, directory2, name, pop, size])

jobs = ['./helperdir/ContamLD_main_DD.R','./helperdir/ContamLD_main_DU.R','./helperdir/ContamLD_main_UU.R','./helperdir/ContamLD_main_full.R']
directory1=sys.argv[1]
directory2=sys.argv[2]
name=sys.argv[3]
pop=sys.argv[4]
size=sys.argv[5]

processes = [Process(target=sample_function, args=([job,directory1,directory2,name,pop,size])) for job in jobs]

for p in processes:
	p.start()

# wait until all processes are finished
for p in processes:
	p.join()
