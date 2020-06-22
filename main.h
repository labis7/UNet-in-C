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
struct maxpoolbackward_data_
{
	/*
	* Conv: The previous result(layer out): the part before we apply maxpooling on it.
	* dpool: 3-dimension image as input of shape: (channels, h, w), its the halfsize(default) of the conv input
	* because its the error that comes backward from smaller to bigger.
	* Output: Size same as conv, we will fill the slots of dconv(output) with the following technique: Find
	* the max element in a sub array of size 2x2 correspoding to the dpool element(error: that came from that 2x2 maxpooling.)
	* and fill this slot with the error, the 3 remaining will be zeros.
	* stride == conv_dimesions == 2
	* Info about current image: channels,dim where dim is the dimension of dpool
	*/
	int channels,dim,o_dim;
	float ***dpool,***conv,***output;
}maxpoolbackward_data;


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

