# myanalysis

To repreduce the memory error: python preselection_memoryerror.py

cpp files of preselection with reduced features added, in order to produce root files including information of first 15 jets. 
To produce the root files(output file path need to be changed, takes several hours): root preselection_forML_reducefeature.cpp 

python scripts to convert the root files with reduced features and information of 15 jets into numpy arrays added.
To produce the npy files(output file path need to be changed): python preselection_root_to_npy.py
