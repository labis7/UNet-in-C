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


struct act_func_data_
{
	int code, channels, dim;
	float ***array;

}act_func_data;


void Activation_Function(struct act_func_data_ *act_func_data)
{
	float ***array;
	int code, channels, dim;
	array = act_func_data->array;
	code = act_func_data->code;
	channels = act_func_data->channels;
	dim = act_func_data->dim;

	float res;
	//Approximation of Sigmoid,  f(x) = x / (1 + abs(x))
	if(code == 1)
	{
		for (int i=0;i<channels;i++)
		{
			for (int j=0;j<dim;j++)
			{
			for (int k=0;k<dim;k++)
				{
					res = (array[i][j][k]/(1 + abs(array[i][j][k])));
					printf("%.2f\t", res);//*(*(*(pA +i) + j) +k));
				}
			printf("\n");
			}
		}
		res = exp(5);
		for (int i=0;i<channels;i++)
		{
			for (int j=0;j<dim;j++)
			{
			for (int k=0;k<dim;k++)
				{
					res = (float)array[i][j][k];
					res = exp(res);
					array[i][j][k] = (1/(1 + exp(- array[i][j][k])));
					printf("%.2f\t", array[i][j][k]);//*(*(*(pA +i) + j) +k));
				}
			printf("\n");
			}
		}
	}
}

float ***make_3darray(int channels,int dim)
{
	int dim1=channels, dim2=dim, dim3=dim;
	int i,j;
	float *** array = (float ***)malloc(dim1*sizeof(float**));

	for (i = 0; i< dim1; i++) {

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
	{
		for (int j=0;j<2;j++)
		{
		for (int k=0;k<2;k++)
			{
				array[i][j][k] = i+j+k;
				//printf("%d\t", array[i][j][k]);//*(*(*(pA +i) + j) +k));
			}
		}
	}


	struct act_func_data_ *ptr_act_func_data = &act_func_data;

	ptr_act_func_data->channels = channels;
	ptr_act_func_data->code = code;
	ptr_act_func_data->dim = dim;
	ptr_act_func_data->array = array;
	Activation_Function(ptr_act_func_data);

	return EXIT_SUCCESS;
}
