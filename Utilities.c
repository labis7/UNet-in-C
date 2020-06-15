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


struct act_func_data_
{
	int code, channels, dim;
	float ***array;

}act_func_data;


void Activation_Function(struct act_func_data_ *act_func_data)
{
	//Aproximation of sigmoid
	//f(x) = x / (1 + abs(x))

}


int main(void) {
	puts("Hello World!"); /* prints Hello World! */


	struct act_func_data_ *ptr_act_func_data = &act_func_data;



	return EXIT_SUCCESS;
}
