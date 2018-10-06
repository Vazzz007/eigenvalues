# invert-reflection
The reflection method for finding the inverse matrix using an adjugate matrix.

Use "make" to compile.

Program supports the following command-line arguments:
  * -i input_file_name.txt - name of the input file
  * -n number - number of elements (default = 10)
  * -v - option for debugging
  * -f formula - define formula (choose from { sym ; symnul ; gilb ; 1 ; 9 }
  * -m number - maximum output size (default = 5)
  
  The samples of using:
  
  make
  
  ./invert -n 200 -f sym -m 7
  
  ./invert -i input.txt
  
  make clean
