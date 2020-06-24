/*
 * Convolutions.c
 *
 *  Created on: Jun 19, 2020
 *      Author: labis
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <main.h>


//////////////////////// CONVOLUTIONS //////////////////////////

void conv(struct conv_data_ *ptr_conv_data)
{
	//f_dim = 3, stride = 1.
	float ***conv_in,***conv_out;
	float ****filter,*bias;
	int dim,f_num,ch_num,o_dim, mode;
	conv_in = ptr_conv_data->conv_in;
	filter = ptr_conv_data->filter;
	bias = ptr_conv_data->bias;
	dim = ptr_conv_data->dim;
	mode = ptr_conv_data->mode;
	ch_num = ptr_conv_data->ch_num;
	f_num = ptr_conv_data->f_num;
	int s=1;
	int f = ptr_conv_data->f_dim; // transp conv :f=2 , conv : f=3

	if (mode >= 1)//padding enabled - probably 1 which means keep the same dim as input
	{
		int pad = mode;
		o_dim = dim;
		int dim_t = dim + 2*pad;
		float ***conv_in_t = make_3darray(ch_num, dim_t);
		conv_out = make_3darray(f_num, o_dim); //number of filters will determine the number of out image channels, dim will be the same in this case.
		//zero padding
		for(int i=0; i< ch_num; i++)
		{
			for(int x = 0; x<pad ; x++)
			{
				for(int y = 0; y< dim_t; y++)
				{
					conv_in_t[i][x][y] = 0;
					conv_in_t[i][y][x] = 0;
					conv_in_t[i][(dim_t-1)-x][y] = 0;
					conv_in_t[i][y][(dim_t-1)-x] = 0;
				}
			}
			//fill the empty center space with conv_in--> then the result wiill be the conv_in padded(conv_in_t)
			for(int x=pad; x<(dim_t-pad); x++)
				for(int y=pad; y<(dim_t-pad); y++)
					conv_in_t[i][x][y] = conv_in[i][x-pad][y-pad];
		}
		// Now we can start the convolution
		int sum;
		for (int i=0; i<f_num; i++)//number of filters
		{
			for(int x=0; x<o_dim; x++)
			{
				for(int y=0; y<o_dim; y++)
				{
					sum=0;
					//seeking on the temp image sub array that we want to mult item wise and then add them for the (x,y) result
					for(int j=0; j < ch_num ; j++)
					{
						for(int k=x; k<(x + f); k++)
						{
							for(int l =y; l<(y+f); l++)
							{
								sum += conv_in_t[j][k][l]*filter[i][j][k-x][l-y];
							}
						}
					}
					conv_out[i][x][y] = sum + bias[i];
				}
			}
		}
	}
	else//mode = pad = 0 = 'normal'
	{
		o_dim = ((dim - f)/s) +1 ;
		conv_out = make_3darray(f_num, o_dim); //number of filters will determine the number of out image channels, dim will be the same in this case.
		int sum;
		for (int i=0; i<f_num; i++)//number of filters
		{
			for(int x=0; x<o_dim; x++)
			{
				for(int y=0; y<o_dim; y++)
				{
					sum=0;
					//seeking on the temp image sub array that we want to mult item wise and then add them for the (x,y) result
					for(int j=0; j < ch_num ; j++)
					{
						for(int k=x; k<(x + f); k++)
						{
							for(int l =y; l<(y+f); l++)
							{
								sum += conv_in[j][k][l]*filter[i][j][k-x][l-y];
							}
						}
					}
					conv_out[i][x][y] = sum + bias[i];
				}
			}
		}

	}

	ptr_conv_data->conv_out = conv_out;
	ptr_conv_data->o_dim = o_dim;
	//number of channels is known before func call,(o_ch)num == f_num)


}




//////////////////// extras ///////////////////////////////////

void crop2half(struct concat_crop_data_ *ptr_concat_crop_data)
{
	float ***image1,***image2,***image3;
	image1 = ptr_concat_crop_data->image1;
	int ch_num = ptr_concat_crop_data->ch_num;
	int dim = ptr_concat_crop_data->dim;
	int o_ch_num = (int)(ch_num/2);

	image2 = make_3darray(o_ch_num, dim);
	image3 = make_3darray(o_ch_num, dim);
	/*
	*starting with the 1st half --> image2,
	* then using the offset o_ch_num+i we
	* can save the other half to the image3 at same time
	*/
	for (int i =0; i < o_ch_num ; i++)
	{
		for (int x =0; x<dim; x++)
		{
			for (int y=0; y<dim; y++)
			{
				image2[i][x][y] = image1[i][x][y];
				image3[i][x][y] = image1[o_ch_num+i][x][y];
			}
		}
	}

	ptr_concat_crop_data->image2 = image2;
	ptr_concat_crop_data->image3 = image3;

	ptr_concat_crop_data->o_ch_num = o_ch_num;

}

void concat(struct concat_crop_data_ *ptr_concat_crop_data)
{
	float ***image1, ***image2, ***image3;
	image1=ptr_concat_crop_data->image1;
	image2=ptr_concat_crop_data->image2;
	int dim = ptr_concat_crop_data->dim;// dimensions for both image1,2 (which is the same)
	int ch_num = ptr_concat_crop_data->ch_num;
	int o_ch_num = ch_num*2;
	image3 = make_3darray(o_ch_num, dim);

	for (int i =0; i < ch_num ; i++)
	{
		for (int x =0; x<dim; x++)
		{
			for (int y=0; y<dim; y++)
			{
				image3[i][x][y] = image1[i][x][y];
				image3[i+ch_num][x][y] = image2[i][x][y];
			}
		}
	}


	ptr_concat_crop_data->image3 = image3;
	ptr_concat_crop_data->o_ch_num = o_ch_num;

}
