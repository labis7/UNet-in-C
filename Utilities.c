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

float ***make_3darray(int channels,int dim);
struct act_func_data_
{
	int code, channels, dim;
	float ***dA, ***Z;

}act_func_data;

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
		for (j = 0; j < dim2; j++) {
			array[i][j] = (float *)malloc(dim3*sizeof(float));
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



	return EXIT_SUCCESS;
}
