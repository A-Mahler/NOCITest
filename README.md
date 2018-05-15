# NOCITest
Test programs for NOCI development

Current files are:  UHFMQC.f03
                    GeneralizedHF.f03
		    NOCI.f03
                    MO.py
                    test_gau.mat
                    
UHFMQC is functional UHF program that will recreate H2 molecular density and energy from .mat file using NOCI formalism
GeneralizedHF is semi-functional GHF program that will recreate density and energy from .mat file of any RHF, UHF, and GHF
  gaussian matrix output
MO.py is utitlity python program for analyzing MO coefficients from rwf dumps

NOCI is primary non-orthogonal CI development code.  Code takes input file and generates psuedo-hamiltonian matrix
  of broken symmetry solutions and diagonalizes matrix to find energy solutions
