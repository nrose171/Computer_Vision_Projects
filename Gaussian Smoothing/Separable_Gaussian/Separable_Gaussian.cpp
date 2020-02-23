#include <iostream>
#include "ReadImage.cpp"
#include "WriteImage.cpp"

#define USE_MATH_DEFINES
#include <cmath>

template<typename T>
void allocateMemory( int x_size, int y_size, T ** matrix );
template<typename T>
void clearMemory( int x_size, T ** matrix );
void Gauss (float s, int Hsize, float * H);
void convolution_1d( bool xHit, int hSize, float * mask, int x_size, int y_size, int ** outputMatrix, int ** inputMatrix );
void innerConvolution_1dx( int hSize, float * mask, int x_size, int x, int y, int ** outputMatrix, int ** inputMatrix );
void innerConvolution_1dy( int hSize, float * mask, int y_size, int x, int y, int ** outputMatrix, int ** inputMatrix );

int main()
{
    int ** input;
    int x_size, y_size, Q;
    char name[20] = "lenna.pgm";
	char outfile_1c[3][20] = {"Sigma=1.pgm", "Sigma=5.pgm", "Sigma=11.pgm"};
	int n = 3;
	int sigma[n] = { 1, 5, 11 };
	int maskScale = 5;

	ReadImage(name, &input, x_size, y_size, Q);
		
	int ** output = new int *[x_size];
 	allocateMemory<int>( x_size, y_size, output );

	for( int i = 0; i<n; i++ )
	{
		int hs = maskScale * sigma[i];

		float * mask = new float[hs];

		Gauss( sigma[i], hs, mask );
		
		convolution_1d( true, hs, mask, x_size, y_size, output, input );
		convolution_1d( false, hs, mask, x_size, y_size, output, output );

		WriteImage( outfile_1c[i], output, x_size, y_size, Q );

		free(mask);
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

//Creates 1D Gaussian Mask
void Gauss (float s, int Hsize, float * H)
{
  
	int     i;
  
	float  cst,tssq,x,sum;
	
	cst=1./(s*sqrt(2.0*M_PI));
  
	tssq=1./(2*s*s);

 	//creation of mask
	for (i=0; i<Hsize; i++) 
	{
    
		x=(float)(i-Hsize/2);

		H[i]=(cst*exp(-(x*x*tssq)));
		//printf( "%lf\n", H[i] );
  
	}
 
  	//printf( "\n" );
	sum=0.0;

	
  	//scaling of mask:
	for (i=0;i<Hsize;i++)
    
		sum += H[i];
  
	for(i=0;i<Hsize;i++)
    {
		H[i] /= sum;
	//printf( "%lf\n", H[i] );
	}

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

//CONVOLUTION FOR ENTIRE 2D INPUT IMG WITH 1D MASK CAN CONVOLVE EITHER X OR Y PIXELS
void convolution_1d( bool xHit, int hSize, float * mask, int x_size, int y_size, int ** outputMatrix, int ** inputMatrix )
{
    for(int row = 0; row < x_size; row++)
	{
        for(int col = 0; col < y_size; col++)
		{
			//convolves with either y pixels or x pixels
            if( xHit == true )
			{
				innerConvolution_1dx( hSize, mask, x_size, row, col, outputMatrix, inputMatrix );
			}
			else if( xHit == false )
			{
				innerConvolution_1dy( hSize, mask, y_size, row, col, outputMatrix, inputMatrix );
			}
        }
    }
}

//INNER CONVOLUTION FOR 2D INPUT IMG WITH 1D MASK FOR X PIXELS
void innerConvolution_1dx( int hSize, float * mask, int x_size, int x, int y, int ** outputMatrix, int ** inputMatrix )
{
	int n = -(hSize / 2);//negative edge for nxn mask
	int p = hSize / 2;//positive edge for nxn mask
	float sum = 0.0;
    for(int i = n; i <= p; i++)
	{
    	if( x-i < 0 || x-i >= x_size )
		{
			sum += 0;
		}
		else
		{
			sum += (mask[i+p]*inputMatrix[x-i][y]);
		}
	}
	if( sum < 0 )
	{
		sum = 0;
	}
	else if( sum > 255 )
	{
		sum = 255;
	}
	outputMatrix[x][y] = sum;
}

//INNER CONVOLUTION FOR 2D INPUT IMG WITH 1D MASK FOR Y PIXELS
void innerConvolution_1dy( int hSize, float * mask, int y_size, int x, int y, int ** outputMatrix, int ** inputMatrix )
{
	int n = -(hSize / 2);//negative edge for nxn mask
	int p = hSize / 2;//positive edge for nxn mask
	float sum = 0.0;
    for(int i = n; i <= p; i++)
	{
    	if( y-i < 0 || y-i >= y_size )
		{
			sum += 0;
		}
		else
		{
			sum += (mask[i+p]*inputMatrix[x][y-i]);
		}
	}
	if( sum < 0 )
	{
		sum = 0;
	}
	else if( sum > 255 )
	{
		sum = 255;
	}
	outputMatrix[x][y] = sum;
}
