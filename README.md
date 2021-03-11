# AC_SingleModule_Tools
Set of tools used for ArgonCube SingleModule run data processing, including light calibration, light reconstruciton, charge/light synchronisation  


SingleModuleLightCalib.C:  
	* Calib(const char* file_light): Update light data file with branches containing number of p.e.  
	
SingleModule3DLightReco.C:  
	* Reconstructed the expected light yield on each light detector using the charge track data.   
	* Output: Fraction of light expected to arrive on each light detector compared to full light yield summed up in 1mm steps.   
	* Output multiplied by scintillation light yield per mm gives total number of photons expected on each light detector   
	* !! NEEDS ROOT 5 !!   
	* run through full datafile:   
		** root -l  
		** .L SingleModule3DLightReco.C  
		** fullrun("my_charge_tracks_file.root")  
		** .q  
	
SingleModule_sync.C:  
	- Synchronises light and charge track events.  
	- Usage:  
		1. root -l  
		2. open(my_light_file, my_charge_file)  
		3. find_offset(), if no peak is visible parameters upper/lower have to be modified  
		4. If a peak is found run find_offset() with n_iter=3 to run algorithm iteratively and automatically call sync function  
		

