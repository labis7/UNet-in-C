/*
 ============================================================================
 Name        : Utilities.c
 Author      : Labiz
 Version     :
 Copyright   : Your copyright notice
 Description : Basic Utilities of a Unet in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>


struct func_data_{
	int (*array)[2][2];
};
struct func_data_ func_data;

/*
 int dim1, dim2, dim3;
 int i,j,k;
 double *** array = (double ***)malloc(dim1*sizeof(double**));

        for (i = 0; i< dim1; i++) {

         array[i] = (double **) malloc(dim2*sizeof(double *));

          for (j = 0; j < dim2; j++) {

              array[i][j] = (double *)malloc(dim3*sizeof(double));
          }

        }

  !!!!!!!!!  OR  !!!!!!!!!!!!!!!!
  int array[dim1][dim2][dim3];
 */
void func(struct func_data_ *p);

int main(void) {
	puts("!!!Hello World!!!"); /* prints !!!Hello World!!! */
	int dim1=2, dim2=2, dim3=2;
	int i,j,k;
	int *** array = (int ***)malloc(dim1*sizeof(int**));

	for (i = 0; i< dim1; i++) {

		array[i] = (int **) malloc(dim2*sizeof(int *));

		for (j = 0; j < dim2; j++) {

			array[i][j] = (int *)malloc(dim3*sizeof(int));
		}

	}
	int x=2,y=2,z=2;
	int (*a)[y][z] = malloc( sizeof(int[x][y][z]) );
	for (int i=0;i<2;i++)
		{
			for (int j=0;j<2;j++)
				{
				for (int k=0;k<2;k++)
					{
						a[i][j][k] = i+k+j+1; //*(*(*(*pt1 +i) + j) +k));//t[i][j][k]);//*(*(*(pA +i) + j) +k));
					}
				}
		}


	int t[2][2][2]= {{{1,2},{3,4}},{{5,6},{7,8}}};;
	int size=2;
	int (*pt)[size][size][size]= &t;
	int (*pt1)[size][size][size];
	pt1=pt;
	//pt = &t;


	for (int i=0;i<2;i++)
	{
		//printf("The elements of A[%d] [j] [k] are given below.\n", i);
		for (int j=0;j<2;j++)
			{
			for (int k=0;k<2;k++)
				{
					//printf("%d\t", (*pt1)[i][j][k]);//*(*(*(*pt1 +i) + j) +k));//t[i][j][k]);//*(*(*(pA +i) + j) +k));
				}
			//printf("\n");
			}
	}
	/* Hint: &t - is of type const pointer to an array of
       2 two dimensional arrays of size [2][2]
       OR
       i can use
       int (*pt)[2][2] = t; //which means t is if type
        const pointer to a 2 dimension array each of them
        consisted of 2d arrays
       */
	//func(pt);
	struct func_data_ *ptr_func_data;
	ptr_func_data = &func_data;
	func_data.array = t;
	func(ptr_func_data);

	return EXIT_SUCCESS;
}

void func(struct func_data_ *p){
	//int t[2][2][2];
	//t= p;
	int (*array)[2][2];
	array = p->array;
	for (int i=0;i<2;i++)
		{
			printf("The elements of A[%d] [j] [k] are given below.\n", i);
			for (int j=0;j<2;j++)
				{
				for (int k=0;k<2;k++)
					{
						printf("%d\t", array[i][j][k]);//*(*(*(pA +i) + j) +k));
					}
				printf("\n");
				}
		}
}
