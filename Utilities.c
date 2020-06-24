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
#include <main.h>

////////////// Function Naming init. ////////////////

////////////////////////////////////////////////////

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

float **make_2darray(int dim1,int dim2)
{
	float **array = (float **)malloc(dim1*sizeof(float*));
	for (int i = 0; i< dim1; i++)
		array[i] = (float *) malloc(dim2*sizeof(float));
	return array;
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

	puts("Main function init!");

	int channels,code,dim;
	code = 1;//sigmoid(approximation)
	channels = 2;
	dim = 5;

	float ***image1;
	image1 = make_3darray(channels,dim);
	for (int i=0;i<channels;i++)
	{
		for (int j=0;j<dim;j++)
		{
			for (int k=0;k<dim;k++)
			{
				image1[i][j][k] = (i+1)*(j*2+k*1) +1;
				printf("%f\t", image1[i][j][k]);//*(*(*(pA +i) + j) +k));
			}
			printf("\n");
		}
		printf("\n");
	}
	printf("\nFILTER:\n");
	float ****filter ;
	int count=1;
	int f_num=4;
	int f=3;
	filter = make_4darray(f_num,channels,f);
	float *bias= (float *)malloc(f_num*sizeof(float));
	for (int i=0;i<f_num;i++)
	{
		for (int l=0; l<channels ; l++)
		{
			for (int j=0;j<f;j++)
			{
				for (int k=0;k<f;k++)
				{
					filter[i][l][j][k] = count;
					count++;
					printf("%f\t", filter[i][l][j][k]);//*(*(*(pA +i) + j) +k));
				}
				printf("\n");
			}
			printf("\n");
		}
		bias[i]=i;
	}


	///////////////////////////////////////////////////////
	////////////////// TESTING SECTION ////////////////////
	/*
	float ***image1;
	image1 = make_3darray(channels,2);
	for (int i=0;i<channels;i++)
		for (int j=0;j<2;j++)
			for (int k=0;k<2;k++)
				image1[i][j][k] = (i+1)*(j*2+k*1);
	*/


	struct conv_data_ *ptr_conv_data = &conv_data;
	ptr_conv_data->conv_in = image1;
	ptr_conv_data->dim = 5;
	ptr_conv_data->mode = 0;
	ptr_conv_data->ch_num=2;
	ptr_conv_data->f_dim=f;
	ptr_conv_data->f_num=f_num;

	ptr_conv_data->bias=bias;
	ptr_conv_data->filter=filter;

	conv(ptr_conv_data);
	float ***res = ptr_conv_data->conv_out;
	int o_dim=ptr_conv_data->o_dim;
	printf("\nResult:\n");
	for (int i=0;i<f_num;i++)
	{
		for (int j=0;j<o_dim;j++)
		{
			for (int k=0;k<o_dim;k++)
			{
				//image1[i][j][k] = (i+1)*(j*2+k*1) +1;
				printf("%f\t", res[i][j][k]);//*(*(*(pA +i) + j) +k));
			}
			printf("\n");
		}
		printf("\n");
	}
	//conv(ptr_conv_data);
	///////////////////////////////////////////////////////
	///////////////////////////////////////////////////////
	printf("\nSuccess! - No Segmentation faults!");
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
