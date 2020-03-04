# JPsiGen
A simple generator that generates exclusive J/psi photo production with it's e-e+ decay mode

********************Istallation***********************
Prerequisite: You need to have a root installed

In the base directory just running "make" should install build the generator `JPsiGen.exe`

You need to add lib directory in the LD_LIBRARY_PATH
For example in .csh enviroment you can do the following command
"setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:`pwd`/lib"

The generator will parse the GenOptions.dat and various parameters of the generator
are set through this file.

One also can use the python script `GenerateGenOptions.py` if you want to 
edit the GenOptions.dat file with command line option, or event run the executable with command line arguments.

An example command is 
`python GenerateGenOptions.py -n 100000 --LUND --Nperfile 35000 -e 10.6 -t -6 --q2Cut 0.02 --ltarg 5 --Run`

* -n 100000             :   Generate 100000 events
* --LUND                :   Write output in the LUND file, otherwise the outpt will be written in the root file
* --Nperfile 35000      :   in a single file number of events should not exceed 35000, as soon 35000 is reach the generator                                 will create andother file, in particular the example command above will create yhree file                                       JPsi_gen_0.txt , JPsi_gen_1.txt and JPsi_gen_2.txt
* -e 10.6               :   The energy of the beam is 10.^ GeV
* -t -6                 :   the upper limit of the Mandelshtam "-t" is 6 GeV^2
* -q2Cut 0.02           :   In the virtual photon flux use only photons with Q2 < 0.02 GeV^2
* -lltarg 5             :   The target length is 5 cm, generated vertexes will be unifirmly distributed in the -2.5 to 2.5                               range
* --Run                 :   Run the generator with specified options
