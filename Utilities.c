/*
 ============================================================================
 Name        : Utilities.c
 Author      : Labiz
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

////////////// Function Naming init. ////////////////
float Random_Normal(int loc, float scale);
float *****testff(float ****t);
float ****make_4darray(int num,int channels,int dim);
float ***make_3darray(int channels,int dim);
////////////////////////////////////////////////////

struct act_func_data_
{
	/*
	 * Z: The output as raw of a layer exactly before we apply the activation function on it.
	 * dA:Thats the difference on the results (after activation function) we get during backprop from the next/forward layer so we
	 * can take advantage of them and backpropagate the error back to activation function backward and then to the previous layer.
	 * Code: Number of commamd: 1)sigmoid, 2)Sigmoid_backward, 3)Relu, 4)Relu_backward
	 * Channels: Number of channel of the specific input matrix
	 * Dim: We assume that we have a square image so the height == width == dim
	 */
	int code, channels, dim;
	float ***dA, ***Z;

}act_func_data;

struct GP_arrays_ //General purpose arrays
{
	float *dim1;
	float **dim2;
	float ***dim3;
	float ****dim4;
	float *****dim5;
}GP_arrays;

struct init_param_
{
	/*
	 * Filters: They made of all the filters for the forward step, including the final 1x1 conv. filters(out_F), These filters 4D dim
	 * as follows: (num_f, num_in_ch, f_h, f_w) and these filters are saved sequencially in a Filters array with the type of(*****)
	 * Bias: Thats the double pointer matric which keeps all the bias values of the network. We need 1-d array for the scalar values
	 * of each bias so the final Bias matrix it is type of (**),so it can include all the different sizes of bias.
	 * F_dc: This matrix contains the decoder upsampling transposed convolution filters.Its type is the same as Filters matrix.
	 */
	int layers, num_f, trim;
	float *****filters,**bias, *****f_dc,**b_dc;

}init_param;

void Initialize_Parameters(struct init_param_ *ptr_init_param)
{
	int layers = ptr_init_param->layers;//number of layers(just the encoder-downsampling number)
	int num_f = ptr_init_param->num_f;  //initial number of filters(they will be doubled for each layer)
	float trim = ptr_init_param->trim;  //weight scale of the values
	float *****filters = (float *****)malloc((layers*2*2-1)*sizeof(float ****));
	float **bias = (float **)malloc((layers*2*2-1)*sizeof(float *));
	float *****f_dc = (float *****)malloc((layers-1)*sizeof(float ****));
	float **b_dc = (float **)malloc((layers-1)*sizeof(float *));

	int ch_in = 1; //number of initial channels input from the input image
	int loc = 0; //normal distribution around zero

	int last_pos=0;
	//Building the downsampling/encoder filters first.(5 layers*2 filters each)
	for (int i=0; i<layers; i++)
	{
		if(i != 0) //after 1st iteration, double the amount of filters of each layer
			num_f = num_f*2; //*16*, 32, 64, 128, 256
		//Building f1
		float ****f1 = make_4darray(num_f, ch_in, 3);//3x3 filter alwars for the simple convolution
		float *b1 = (float *)malloc(num_f*sizeof(float));
		for (int x=0;x<num_f;x++)
		{
			for (int y=0;y<ch_in;y++)
				for (int z=0;z<3;z++)
					for(int w=0;w<3;w++)
					{
						f1[x][y][z][w] = Random_Normal(loc, trim);
					}
			b1[x] =  Random_Normal(loc, trim);
		}

		/*
		 * Very important! Next filter will have as input the output/channels of the previous filter
		 * More : When we apply a number of filters on an image, we create num_f channels on the image result.
		 * So each filter must have input dim same as the input image channels, output dim = channels we want to occur after conv.
		 */
		ch_in = num_f;

		float ****f2 = make_4darray(num_f, ch_in, 3);//3x3 filter alwars for the simple convolution
		float *b2 = (float *)malloc(num_f*sizeof(float));
		for (int x=0;x<num_f;x++)
		{
			for (int y=0;y<ch_in;y++)
				for (int z=0;z<3;z++)
					for(int w=0;w<3;w++)
					{
						f2[x][y][z][w] = Random_Normal(loc, trim);
					}
			b2[x] = Random_Normal(loc, trim);
		}
		filters[2*i] = f1;
		filters[2*i +1] = f2;
		bias[2*i] = b1;
		bias[2*i +1] = b2;
		last_pos = i; //last position of the i that shows the current layer, we will need it later so we keep saving filters and
		//bias sequencially.
	}
	for (int i = 1 ; i < layers; i++)
	{
		num_f = (int)num_f/2;//It will be always power of number of filters,so the division will give back integer

		float ****fdc = make_4darray(num_f, ch_in, 2);//2x2 filter always for the upsampling/transposed convolution
		float *bdc = (float *)malloc(num_f*sizeof(float));
		float ****f1 = make_4darray(num_f, ch_in, 3);//3x3 filter always for the simple convolution
		float *b1 = (float *)malloc(num_f*sizeof(float));
		for (int x=0;x<num_f;x++)
		{
			for (int y=0;y<ch_in;y++)
			{
				for (int z=0;z<3;z++)
				{
					for(int w=0;w<3;w++)
						f1[x][y][z][w] = Random_Normal(loc, trim);
				}
				for (int z=0;z<2;z++)
				{
					for(int w=0;w<2;w++)
						fdc[x][y][z][w] = Random_Normal(loc, trim);
				}
			}
			bdc[x] =  Random_Normal(loc, trim);
			b1[x] = Random_Normal(loc, trim);
		}

		ch_in = num_f ;

		float ****f2 = make_4darray(num_f, ch_in, 3);//3x3 filter always for the simple convolution
		float *b2 = (float *)malloc(num_f*sizeof(float));
		for (int x=0;x<num_f;x++)
		{
			for (int y=0;y<ch_in;y++)
				for (int z=0;z<3;z++)
					for(int w=0;w<3;w++)
					{
						f2[x][y][z][w] = Random_Normal(loc, trim);
					}
			b2[x] =  Random_Normal(loc, trim);
		}

		filters[last_pos*2] = f1;
		filters[last_pos*2 + 1] = f2;
		bias[last_pos*2] = b1;
		bias[last_pos*2 + 1] = b2;
		f_dc[i-1] = fdc;
		b_dc[i] = bdc;
		last_pos++;
	}
	ptr_init_param->filters = filters;
	ptr_init_param->bias = bias;
	ptr_init_param->f_dc = f_dc;
	ptr_init_param->b_dc = b_dc;
}

float ***Activation_Function(struct act_func_data_ *act_func_data)
{
	float ***Z;
	int code, channels, dim;
	Z = act_func_data->Z;
	code = act_func_data->code;
	channels = act_func_data->channels;
	dim = act_func_data->dim;

	float ***res = make_3darray(channels, dim);
	//Approximation of Sigmoid,  f(x) = x / (1 + abs(x))
	if(code == 1)
	{
		for (int i=0; i<channels; i++)
			for (int j=0; j<dim; j++)
				for (int k=0; k<dim; k++)
					res[i][j][k] = (1/(1 + exp(- Z[i][j][k]))); //Thats the original sigmoid!
	}
	else if(code == 2)//Sigmoid backpropagation function
	{
		float ***dA;
		dA =act_func_data->dA;
		for (int i=0; i<channels; i++)
			for (int j=0; j<dim; j++)
				for (int k=0; k<dim; k++)
					res[i][j][k] = (1/(1 + exp(- Z[i][j][k])));//find sig
		for (int i=0; i<channels; i++)
			for (int j=0; j<dim; j++)
				for (int k=0; k<dim; k++)
					res[i][j][k] = (dA[i][j][k])*(res[i][j][k])*(1 - res[i][j][k]);
	}
	else if(code == 3) //RELU activation function
	{
		for (int i=0; i<channels; i++)
			for (int j=0; j<dim; j++)
				for (int k=0; k<dim; k++)
				{
					if(Z[i][j][k] <= 0)
						res[i][j][k] = 0;
					else
						res[i][j][k] = Z[i][j][k];
				}
	}
	else //Relu backpropagation function
	{
		float ***dA;
		dA = act_func_data->dA;
		for (int i=0; i<channels; i++)
			for (int j=0; j<dim; j++)
				for (int k=0; k<dim; k++)
				{
					if(Z[i][j][k] <= 0)
						res[i][j][k] = 0;
					else
						res[i][j][k] = dA[i][j][k];
				}
	}
	return res;
}

float ***make_3darray(int channels,int dim)
{
	int dim1=channels, dim2=dim, dim3=dim;
	int i,j;
	float *** array = (float ***)malloc(dim1*sizeof(float**));

	for (i = 0; i< dim1; i++)
	{
		array[i] = (float **) malloc(dim2*sizeof(float *));
		for (j = 0; j < dim2; j++)
			array[i][j] = (float *)malloc(dim3*sizeof(float));
	}
	return array;
}
float ****make_4darray(int num,int channels,int dim)
{
	int dim0=num, dim1=channels, dim2=dim, dim3=dim;
	int i,j,k;
	float **** array = (float ****)malloc(dim0*sizeof(float***));

	for (i = 0; i< dim0; i++)
	{
		array[i] = (float ***) malloc(dim1*sizeof(float **));
		for (j = 0; j < dim1; j++) {
			array[i][j] = (float **)malloc(dim2*sizeof(float *));
			for (k = 0; k < dim2; k++)
				array[i][j][k] = (float *)malloc(dim3*sizeof(float));
		}
	}
	return array;
}

float Random_Normal(int loc, float scale)
{
   // return a normally distributed random value
	scale = 1;
	loc=0;
	float v1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. ); //random gen
	float v2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );//random gen
	return ((cos(2*3.14*v2)*sqrt(-2.*log(v1)))*scale + loc);
}
int main(void) {
	time_t t;
	srand((unsigned) time(&t));

	puts("Hello World!"); /* prints Hello World! */

	int channels,code,dim;
	code = 1;//sigmoid(approximation)
	channels = 2;
	dim = 2;

	float ***array;
	array = make_3darray(channels,dim);
	for (int i=0;i<2;i++)
		for (int j=0;j<2;j++)
			for (int k=0;k<2;k++)
				array[i][j][k] = i+j+k;
				//printf("%d\t", array[i][j][k]);//*(*(*(pA +i) + j) +k));


	struct init_param_ *ptr_init_param = &init_param;

	//Initialize_parameters(ptr_init_param);
	ptr_init_param->layers = 5;
	ptr_init_param->num_f = 16;
	ptr_init_param->trim = 0.1;
	printf("yo\n");
	Initialize_Parameters(ptr_init_param);
	float *****filters = ptr_init_param->filters;
	float *****f_dc = ptr_init_param->f_dc;
	float **bias= ptr_init_param->bias;
	float **b_dc=ptr_init_param->b_dc;


	int max_layers = 10;
	int filter_layer = 1;
	float ****f_t1,****f_t2;
	f_t1 = filters[filter_layer - 1];
	f_t2 = filters[filter_layer];
	// Finding dimensions of the layer filters using only layer number
	int num_filt,input_ch,f_dim;
	num_filt = 16;
	input_ch =1;
	f_dim=3;
	/////////////////////////
	for (int i=0;i<num_filt;i++)
	{
		for (int j=0;j<input_ch;j++)
		{
			for (int k=0;k<f_dim;k++)
			{
				printf("\n");
				for (int l=0;l<f_dim;l++)
				{
					printf("%.2f\t",f_t1[i][j][k][l]);
				}
			}
		}
	}
	return EXIT_SUCCESS;
}

float *****testff(float ****t){
	float ****t2 = make_4darray(1,2,3);
	for (int i=0;i<1;i++)
		for (int j=0;j<2;j++)
			for (int k=0;k<3;k++){
				for(int y=0;y<3;y++)
				{
					t[i][j][k][y]=i+j+y+k;
					t2[i][j][k][y]=2*i+2*j+2*y+2*k;
				}
			}
	float *****t_array = (float *****)malloc(2*sizeof(float ****));
	t_array[0] = t;
	t_array[1] = t2;
	return t_array;
}
