# pw2cif
Convert Quantum Espresso VC-relax output to a CIF file. 

Process Quantum Espresso pw.x output file into a CIF file. 

Ported from perl script qe2cif.pl by Alexandr Fonari. 
Updated to read atomic positions from current version of QE pw.x output file. 

Author: mikas.remeika@bk.tsukuba.ac.jp
 
Tested with 'vc-relax', QE v6.2.
 
Input parameters:
 
Which atomic positions to process: -sf (start and final) -int (intermediate).
Intermediate positions do not include start and final positions.
All positions will be written into one CIF file, labeles accordingly. 

pw.x output file (captuted output to stdout). 
