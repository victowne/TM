import numpy as np 
import h5py

#create hdf5 file
f = h5py.File('groups.hdf5','w-') #r for read only r+ for read and write 
#create group,dataset and attributes
subgroup = f.create_group('subgroup')
subsubgroup = subgroup.create_group('ssubgruop')
f.create_group('grp1/grp2')

f['data1'] = 5.0
f['data2'] = [1,2,3,4,5,0,6,5]
f['data3'] = 3
f.create_dataset('grp1/grp2/dset1',data=1.0)

dset = f['grp1/grp2/dset1'] #only dataset have attributes
dset.attrs['value'] = 2134652.5
dset.attrs['title'] = 'kfc'

###### visit #######
def printname(name):
    print name
f.visit(printname)
####################