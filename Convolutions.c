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
