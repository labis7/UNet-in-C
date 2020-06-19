/*
 ============================================================================
 Name        : Utilities.c
 Author      : Labiz
 Version     : V0.5
 Copyright   : Your copyright notice
 Description : Utilities of U-Net, including initializations and important 3d,4d matrix building tools
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
void hello();
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
struct init_GN_
{
	/*
	 *-Layers: We need the number of forward layers(like in init_param) so we can calculate the final number of elements
	 *gamma/beta matrix will have.
	 *-Starting_num_ch: This is the number of filters we apply in the very first convolution of the net, which is always 16,
	 *that way we will able to calculate the rest dimension of gamma/beta sub-matrices.
	 *-Gamma/Beta: These two 'big' matrices are going to keep all the elements(sub-arrays) of each layer, so we can access them later.
	 */
	int layers, starting_num_ch, trim;
	float **gamma, **beta;

}init_GN;

void Initialize_GN(struct init_GN_ *ptr_init_GN)//Group Normalization Init.
{
	int layers, num_ch, trim;
	layers = ptr_init_GN->layers; //We PUT THE FORWARD number of layers!! so it can match the other ini func
	trim = ptr_init_GN->trim; //optimal setting:0.05.
	num_ch = ptr_init_GN->starting_num_ch;//DEFAULT:16
	float **gamma = (float **)malloc((layers*2*2-1)*sizeof(float *));
	float **beta = (float **)malloc((layers*2*2-1)*sizeof(float *));

	//Following the same logic as the init_param function:
	//First we build the forward pass gamma/beta, then the upsampling and finally the out gamma/beta(10th layer).
	int last_pos=0,loc=0;
	for (int i=0; i<layers; i++)
	{
		if(i != 0) //after 1st iteration, double the amount of filters of each layer
			num_ch = num_ch*2; //*16*, 32, 64, 128, 256
		float *gamma1 = (float *)malloc(num_ch*sizeof(float));
		float *gamma2 = (float *)malloc(num_ch*sizeof(float));
		float *beta1 = (float *)malloc(num_ch*sizeof(float));
		float *beta2 = (float *)malloc(num_ch*sizeof(float));
		for (int x=0; x < num_ch; x++)
		{
			gamma1[x] =  Random_Normal(loc, trim);
			gamma2[x] =  Random_Normal(loc, trim);
			beta1[x] =  Random_Normal(loc, trim);
			beta2[x] =  Random_Normal(loc, trim);
		}

		gamma[2*i] = gamma1;
		gamma[2*i+1] = gamma2;
		beta[2*i] = beta1;
		beta[2*i+1] = beta2;
		last_pos++;
	}

	for (int i=1; i<layers; i++)
	{
		num_ch = (int)num_ch/2;//It will be always power of number of filters,so the division will give back integer

		float *gamma1 = (float *)malloc(num_ch*sizeof(float));
		float *gamma2 = (float *)malloc(num_ch*sizeof(float));
		float *beta1 = (float *)malloc(num_ch*sizeof(float));
		float *beta2 = (float *)malloc(num_ch*sizeof(float));
		for (int x=0; x < num_ch; x++)
		{
			gamma1[x] =  Random_Normal(loc, trim);
			gamma2[x] =  Random_Normal(loc, trim);
			beta1[x] =  Random_Normal(loc, trim);
			beta2[x] =  Random_Normal(loc, trim);
		}

		gamma[2*last_pos] = gamma1;
		gamma[2*last_pos+1] = gamma2;
		beta[2*last_pos] = beta1;
		beta[2*last_pos+1] = beta2;
		last_pos++;
	}
	//Now we build last layer/output group normalization paramenters
	num_ch = 1;

	float gamma1; //= (float *)malloc(num_ch*sizeof(float));
	float beta1; //= (float *)malloc(num_ch*sizeof(float));

	gamma1 = Random_Normal(loc, trim);
	beta1 =	Random_Normal(loc, trim);

	gamma[last_pos*2] = &gamma1;
	beta[last_pos*2] = &beta1;

	//Filling up the 'big' arrays of all the layer data
	ptr_init_GN->gamma = gamma;
	ptr_init_GN->beta = beta;


}


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
		last_pos++; //last position of the i that shows the current layer, we will need it later so we keep saving filters and
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
	printf("\nlast_pos:%d",last_pos);
	//Now build the last one 1x1 conv filters
	num_f = 1; //we need 1 channels for the output(same as input)
	float ****out_f = make_4darray(num_f, ch_in, 1);//1x1 output filter
	float *out_b = (float *)malloc(num_f*sizeof(float));
	for (int x=0;x<num_f;x++)
	{
		for (int y=0;y<ch_in;y++)
			for (int z=0;z<1;z++)
				for(int w=0;w<1;w++)
				{
					out_f[x][y][z][w] = Random_Normal(loc, trim);
				}
		out_b[x] =  Random_Normal(loc, trim);
	}
	filters[last_pos*2] = out_f;
	bias[last_pos*2] = out_b;

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
	///////////////////////////////////////////////////////
	////////////////// TESTING SECTION ////////////////////



	///////////////////////////////////////////////////////
	///////////////////////////////////////////////////////
	printf("\nSuccess!");
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
