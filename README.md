# eigenvalues_QR_reflection
The QR reflection method for finding the eigenvalues using an QR decomposition.

Use "make" to compile.

Program supports the following command-line arguments:
  * -i input_file_name.txt - name of the input file
  * -n number - number of elements (default = 10)
  * -v - option for debugging
  * -f formula - define formula (choose from { 10 } // from the paper
  * -m number - maximum output size (default = 5)
  
  The samples of using:
  
  make
  
  ./invert -n 1000 -f 10 -m 10
  
  ./invert -i input.txt -v
  
  for i in 2 3 4 5 6 7 8 9; do ./invert -i input2_$i.txt; done
  
  make clean
