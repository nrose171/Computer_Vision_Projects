#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include "ReadImage.cpp"
#include "WriteImage.cpp"

using namespace std;

void Gauss (float s, int Hsize, float * H);
void convolution_1d_image( int hSize, float * mask, int x_size, float * outputMatrix, float * inputMatrix );
void inner_convolution_1d_image( int hSize, float * mask, int x_size, int x, float * outputMatrix, float * inputMatrix );
void read_txt_file( char fName[], float * inputMatrix, int & x_size );
void write_txt_file( char fName[], float * inputMatrix, int x_size );

int main()
{
    float * input = new float [500]; //Max txt file size
    int x_size, y_size, Q;
    char name[20] = "Rect_128.txt";
    char outfile[3][20] = {"1d_sigma=1.txt", "1d_sigma=5.txt", "1d_sigma=11.txt"} ; //names of output files
	int n = 3;	//Number of times to smooth input file
	int sigma[n] = {1, 5, 11};	//Set sigma values for mask here

	read_txt_file( name, input, x_size );

	for( int i = 0; i < n; i++ )
	{
		int hs = 5*sigma[i];
		float * mask = new float [hs];
		Gauss( sigma[i], hs, mask );

		float * output = new float [x_size];
		convolution_1d_image( hs, mask, x_size, output, input );
		write_txt_file( outfile[i], output, x_size );

		free(mask);
		free(output);
	}

	free(input);

    return 0;
}

void read_txt_file( char fName[], float * inputMatrix, int & x_size )
{
	ifstream inputFile;
	int counter = 0;

	inputFile.open(fName);
	string strArr[500];
	string::size_type sz;
	
	string line;

	if( inputFile.is_open() )
	{
		while( getline( inputFile, line ) )
		{
			istringstream iss(line);
			if( (iss >> strArr[counter]) )
			{
				//cout << counter << endl;
				inputMatrix[counter] = stof( strArr[counter], &sz );
				//cout << input[counter] << endl;
				counter++;
				continue;
			}
		}
	}
	else
	{
		cout << "Unable to open" << endl;
	}
	//cout << "finished" << endl;
	inputFile.close();
	x_size = counter;
}

void write_txt_file( char fName[], float * inputMatrix, int x_size )
{
	ofstream of;
	of.open( fName );

	for( int i = 0; i < x_size; i++ )
	{
		of << inputMatrix[i] << endl;
	}
	of.close();
}

//CONVOLUTION FOR ENTIRE 1D INPUT IMG WITH 1D MASK
void convolution_1d_image( int hSize, float * mask, int x_size, float * outputMatrix, float * inputMatrix )
{
	for( int x = 0; x < x_size; x++ )
	{
		inner_convolution_1d_image( hSize, mask, x_size, x, outputMatrix, inputMatrix );
	}
}

//INNER CONVOLUTION FOR 1D INPUT IMG AT PIXEL(X) WITH 1D MASK
void inner_convolution_1d_image( int hSize, float * mask, int x_size, int x, float * outputMatrix, float * inputMatrix )
{
	int p = hSize / 2;
	int n = -p;
	float sum = 0.0;

	for( int i = n; i <= p; i++ )
	{
		//if mask convolution is outside of bounds, replace sum with 0
		if( x-i < 0 || x-i >= x_size )
		{
			sum += 0;
		}
		else
		{
			sum += (mask[i+p]*inputMatrix[x-i]);
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
	outputMatrix[x] = sum;
}

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
