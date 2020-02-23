#include <iostream>
#include <string>
#include <cmath>
#include "ReadImage.cpp"
#include "WriteImage.cpp"

template<typename T>
void allocateMemory( int x_size, int y_size, T ** matrix );
template <typename T>
void clearMemory( int x_size, T ** matrix );
void gaussMask_2d( int hSize, float sigma, float ** mask );
void convolution_2d( int hSize, float ** mask, int x_size, int y_size, int ** outputMatrix, int ** inputMatrix );
void innerConvolution_2d( int hSize, float ** mask, int x_size, int y_size, int x, int y, int ** outputMatrix, int ** inputMatrix );

int main()
{
    int ** input;
    int x_size, y_size, Q;
    char name[20] = "lenna.pgm";	//input file name
    char outfile_2d[3][20] = {"2D_Sigma=1.pgm", "2D_Sigma=5.pgm", "2D_Sigma=11.pgm"};	//output file names
	int n = 3;	//Number of smoothings done
	int sigma[n] = { 1, 5, 11 };	//Different sigma levels used
	int maskScale = 5;

	ReadImage(name, &input, x_size, y_size, Q);
		
	int ** output = new int *[x_size];
 	allocateMemory<int>( x_size, y_size, output );
	
	for( int i = 0; i < n; i++ )
	{
		int hs = maskScale * sigma[i];
		float ** mask = new float *[hs];
		allocateMemory<float>( hs, hs, mask );

		gaussMask_2d( hs, sigma[i], mask );
		convolution_2d( hs, mask, x_size, y_size, output, input );
		WriteImage(outfile_2d[i], output, x_size, y_size, Q);

		clearMemory<float>( hs, mask );
	}

	clearMemory<int>( x_size, output );

    return 0;
}

//ALLOCATES MEMORY FOR A PASSED 2D MATRIX
template<typename T>
void allocateMemory( int x_size, int y_size, T ** matrix )
{
	for( int i = 0; i < x_size; i++ )
	{
		matrix[i] = new T[y_size];
	}
}

//FREES MEMORY FOR A PASSED 2D MATRIX
template <typename T>
void clearMemory( int x_size, T ** matrix )
{
	for( int i = 0; i < x_size; i++ )
	{
		free( matrix[i] );
	}
	free( matrix );
}

//CREATES A SQUARE 2D GAUSSIAN MASK WITH SUM OF 1 SO THAT EVERY ELEMENT IS A PROBABILITY
void gaussMask_2d( int hSize, float sigma, float ** mask )
{
	
	float x, y;
	
	float cst = 1/(2*M_PI*sigma*sigma);
	float tssq = 1/(2*sigma*sigma);
	
	float sum = 0;
	for( int rows = 0; rows < hSize; rows++ )
	{
		x=(float)(rows-hSize/2);
		for( int cols = 0; cols < hSize; cols++ )
		{
			y = (float)(cols-hSize/2);
			mask[rows][cols] = cst * exp( -(((x*x) + (y*y)) * tssq) );
			sum += mask[rows][cols];
		}
	}

	float newSum = 0;
	for( int rows = 0; rows < hSize; rows++ )
	{
		for( int cols = 0; cols < hSize; cols++ )
		{
			mask[rows][cols] /= sum;
			newSum += mask[rows][cols];
		}
	}
}

//CONVOLUTION FOR ENTIRE 2D INPUT IMG WITH 2D MASK
void convolution_2d( int hSize, float ** mask, int x_size, int y_size, int ** outputMatrix, int ** inputMatrix )
{
	//CALLS INNER CONVOLUTION FOR EACH PIXEL
	for( int y = 0; y < y_size; y++ )
	{
		for( int x = 0; x < x_size; x++ )
		{
			innerConvolution_2d( hSize, mask, x_size, y_size, x, y, outputMatrix, inputMatrix );
		}
	}
}

//INNER CONVOLUTION FOR 2D INPUT IMG WITH 2D MASK AT PIXEL(X,Y) 
void innerConvolution_2d( int hSize, float ** mask, int x_size, int y_size, int x, int y, int ** outputMatrix, int ** inputMatrix )
{
	int n = -(hSize / 2);//negative edge for nxn mask
	int p = hSize / 2;//positive edge for nxn mask

	int hCheck = hSize % 2;
	if( hCheck == 0 )
	{ 
		float sum = 0;
		for( int k = n; k < p; k++ )
		{
			for( int j = n; j < p; j++ )
			{
				//if mask convolution is outside of bounds, replace sum with 0
				if( ((y-j) < 0) || ((x-k) < 0) || ((y-j) >= y_size) || ((x-k) >= x_size) )
				{
					sum += 0;
				}
				else
				{
					sum += (mask[j+p][k+p] * inputMatrix[y-j][x-k]);
				}	
			}
		}
		//Normalize pixels
		if( sum < 0 )
		{
			sum = 0;
		}
		else if( sum > 255 )
		{
			sum = 255;
		}
		outputMatrix[y][x] = sum;
	}
	else
	{
		float sum = 0;
		for( int k = n; k <= p; k++ )
		{
			for( int j = n; j <= p; j++ )
			{
				//if mask convolution is outside of bounds, replace sum with 0
				if( ((y-j) < 0) || ((x-k) < 0) || ((y-j) >= y_size) || ((x-k) >= x_size) )
				{
					sum += 0;
				}
				else
				{
					sum += (mask[j+p][k+p] * inputMatrix[y-j][x-k]);
				}
				
			}
		}
		//Normalize pixels
		if( sum < 0 )
		{
			sum = 0;
		}
		else if( sum > 255 )
		{
			sum = 255;
		}
		outputMatrix[y][x] = sum;
	}
	

}

