/*
 * main.h
 *
 *  Created on: Jun 20, 2020
 *      Author: labis
 */

#ifndef MAIN_H_
#define MAIN_H_

float Random_Normal(int loc, float scale);
float *****testff(float ****t);
float ****make_4darray(int num,int channels,int dim);
float ***make_3darray(int channels,int dim);

#endif /* MAIN_H_ */
struct maxpool_data_
{
	/*
	* Image: 3-dimension image as input of shape: (channels, h, w)
	* the output size will be OH = (int)(H - 2)/2 + 1,
	* it will be always integer result since we use resolution in the power of 2.
	* stride == conv_dimesions == 2
	* Info about current image: channels,dim
	*/
	int channels,dim,o_dim;
	float ***image,***output;
}maxpool_data;
