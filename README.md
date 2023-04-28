The GitHub repository contains a C++ tool for estimating upper bounds on code (Hamming) distance estimation for QC codes using an MPICH parallel implementation. The tool supports multithreaded and server clusters and is based on the general idea proposed by MacKay, which has been further developed by Smarandache-Vontobel and improved by Butler-Siegel.

This implementation uses the NTL library, MPI CH, and includes compiled libraries and MS Visual Studio projects. The tool provides a valuable resource for researchers and practitioners interested in optimizing the performance of their error-correcting coding schemes.

Estimating the code distance is a crucial step in designing efficient error-correcting codes that can recover from errors introduced during communication. By providing an estimation of the upper bound on the code distance, this tool allows researchers and practitioners to optimize their coding schemes for better error-correction performance.

Overall, this repository offers a comprehensive C++ tool for estimating upper bounds on code distance estimation for QC codes. With its support for multithreading and server clusters, researchers and practitioners will find this tool valuable in exploring various strategies for designing efficient error-correcting codes with high-performance rates.





Example of use in bin folder (for 32 treads): mpiexec.exe -n 32 MSVBS_bound.exe -inputFile code_to_check.txt -outputFile out.txt

code_to_check.txt contain list of file's names with QC-LDPC  parity-check matrix.

in bin bat file contain example of code distance estimation for (3,6) regular QC-LDPC code and  WIFI 802.11ad (3x16x42, 4x16x42 , 8x16x42, 6x16x42) QC-LDPC codes.

Result of estimation out.txt:
Code distance name_of_file 

24 3x6x30.txt


6 3x16x42.txt


10 4x16x42.txt


17 8x16x42.txt


14 6x16x42.txt

Format of file
Variable Nodes numbers| Check Nodes numbers in protograph | CPMs size
after circulant shifts list
3x6x30.txt:


6 3 30


8 9 28 7 5 7	

25 8 6 23 3 0


2 12 5 26 9 15	

8x16x42.txt:


16 8 42


40 -1 38 -1 13 -1 5 -1 18 -1 -1 -1 -1 -1 -1 -1


34 -1 35 -1 27 -1 -1 30 2 1 -1 -1 -1 -1 -1 -1


-1 36 -1 31 -1 7 -1 34 -1 10 41 -1 -1 -1 -1 -1


-1 27 -1 18 -1 12 20 -1 -1 -1 15 6 -1 -1 -1 -1


35 -1 41 -1 40 -1 39 -1 28 -1 -1 3 28 -1 -1 -1


29 -1 0 -1 -1 22 -1 4 -1 28 -1 27 -1 23 -1 -1


-1 31 -1 23 -1 21 -1 20 -1 -1 12 -1 -1 0 13 -1


-1 22 -1 34 31 -1 14 -1 4 -1 -1 -1 13 -1 22 24


Reference
1. D. J. MacKay and M. C. Davey, “Evaluation of Gallager codes for short block length and high rate applications,” Proc. of the IMA Workshop on Codes, System and Graphical Models, 1999. Springer-Verlag 2001, pp. 113–130
2. R. Smarandache and P. O. Vontobel, “Quasi-cyclic LDPC codes: Influence of proto- and Tanner-graph structure on minimum Hamming distance upper bounds,” IEEE Trans. Inf. Theory, vol. 58, no. 2, pp. 585–607, Feb. 2012
3. Brian K. Butler, Paul H. Siegel Bounds on the Minimum Distance of Punctured Quasi-Cyclic LDPC Codes. IEEE Transactions on Information Theory 2013

Source codes for fast estimation of permanent (c and matlab) could be taked from Brian K. Butler site: https://sites.google.com/site/brianbutlerengineer/home/research
