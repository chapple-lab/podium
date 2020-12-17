# PODIUM
Modified version of XCMS to handle SIL detection for mapping entire metabolic networks
This package is an extension of xcms and includes a new object, xcmsSet2.

This object holds all xcmsRaw files in memory (unlike xcmsSet) which can speed computation by eliminating
the need to constantlty read in xcmsRaw files whenever they are called by a method.

The xcmsSet2 object also includes some modified methods, most notably diffReport and diffReport2, that have
been adapted from base xcms to accomodate the new information generated by our experimental approach.
The majority of the methods responsible for the SIL analysis (eg. the pairing algorithm), however, are not included in this package.
They have been instead grouped in the [podiumMethods](https://github.com/chapple-lab/podiumMethods) repository for the time being.

# Under Construction
This repository is currently under construction.  We hope to have streamlined documentation and a nicer landing page available soon.  For now, you can find installation and setup instructions [in the supplemental document](https://github.com/chapple-lab/podium/blob/master/buildignore/Supplement_and_Instruction_Manual.pdf).  You can download the latest release of PODIUM [here](https://github.com/chapple-lab/podium/releases/download/initial_pre_release/podium_1.2.0.tar.gz) and PODIUMmethods [here](https://github.com/chapple-lab/podiumMethods/raw/master/buildignore/sources/podiumMethods_1.2.0.tar.gz).
