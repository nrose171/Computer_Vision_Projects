
Basic Gaussian Smoothing Projects:

Need g++ compiler installed.

1D_Gaussian)

	Input needed:
	.txt file (Ex. Rect_128.txt)

	To Compile:
	g++ -o <executable filename> 1D_Gaussian.cpp
	
	To Run:
	./<executable filename>

	Output:
	n .txt files that show a gaussian blur with n seperate sigma values applied to the input file.

2D_Gaussian)
	
	Input needed:
	.pgm file (Ex. lenna.pgm)
	
	To Compile:
	g++ -o <executable filename> 2D_Gaussian.cpp

	To Run:
	./<executable filename>

	Output:
	n .pgm files with gaussian blurs of different sigmas applied to the input file.

Seperable_Gaussian)

	Input:
	.pgm file (Ex. lenna.pgm)

	To Compile:
	g++ -o <executable filename> separable_Gaussian.cpp

	To Run:
	./<executable filename>

	Outputs:
	n .pgm files with gaussian blurs of different sigmas applied to the input file.
	
