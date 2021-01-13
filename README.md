# inner-shell
laser rescattering excitation of atomic inner shell

__HOW TO RUN THIS WHOLE PROCESS__

!(campus_map2.jpg)

1. Select an element for this process

2. Using rksuite.cpp, rksuite.h, and IonProbLi.cpp (changes will have to be made for atomic parameters of the selected element), generate the ion populations for the
selected element, all the way down to the bare nucleus. If running on the Caviness cluster, you will also need the job file jobLi.qs. Once these curves are generated, plot
them up and determine which ion number will be your starting state (an ion with an empty outer shell, roughly 90% population) and which ion number will be your final ionization
state (the next ionized species, at roughly 10% ionization.) This code should output a .dat file.

3. Using TableGammaLambdaI.nb, find the appropriate wavelength and intensity combinations (will be generated by functions in the file once specific parameters are inputted.)

4. Using rksuite.cpp, rksuite.h, ran.h, nr3.h, deviates.h, and cutoff_v5.cpp, run the cutoff code as the intermediary step. This code will require some atomic parameters as well
as the intensities generated in the previous step. This code should output a .dat file. If running on the Caviness cluster, you will also need the job file jobCutoff.qs.

5. Taking the .dat file from the cutoff code and using the files logger.h, logger.cpp, data1D.h, data1D.cpp, collection1D.h, collection1D.cpp, and postProcess_v3.cpp, run the
post processing of the data generated by the cutoff process. There will be a few things to fill in in the postProcess_v3.cpp file about the naming of the input files generated
previously. If running on the cluster, you will need jobPostProcess.qs. This code will output two files; the one containing the useful flux data is the binNorm file.

6. With the binNorm file and the cross-section data for the element in question (NIST or other resources), the fluence can be calculated by multiplying the flux and the
cross-section and integrating.
