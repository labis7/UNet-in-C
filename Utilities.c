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
float *****testff(float ****t);
void Initialize_parameters(struct init_param_ *ptr_init_param);
float ***make_3darray(int channels,int dim);
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
	float *****filters,**bias, *****f_dc;

}init_param;

void Initialize_parameters(struct init_param_ *ptr_init_param)
{

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

int main(void) {
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


	struct act_func_data_ *ptr_act_func_data = &act_func_data;

	ptr_act_func_data->channels = channels;
	ptr_act_func_data->code = code;
	ptr_act_func_data->dim = dim;
	ptr_act_func_data->Z = array;
	float ***res = Activation_Function(ptr_act_func_data);

	float ****test = make_4darray(1,2,3);
	//struct act_func_data_ *ptr_init_params = &init_params;

	return EXIT_SUCCESS;
}

