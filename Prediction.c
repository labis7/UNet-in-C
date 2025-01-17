#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <main.h>


void predict(struct images_data_ *images_data,struct params_ *params)
{

	printf("Choose an Image(number 0-%d) for prediction: ",(images_data->im_num -1));
	int predict_num;
	scanf("%d", &predict_num);

	printf("\n");

	//choose image for prediction(it can be a list saved in the struct)
	float ****images=images_data->images;
	float ***image = images[predict_num];//1st image is selected
	float ****labels=images_data->labels;
	float ***label = labels[predict_num];
	int dim = images_data->dim; // same across the whole data
	///////// Laod Parameters  //////////
	float *****filters = params->filters;
	float **bias = params->bias;
	float *****f_dc = params->f_dc;
	float **b_dc = params->b_dc;
	//int gn_batch = params->gn_batch;
	//float **gamma = params->ga;
	//float **beta = params->be;
	//int layers =params->layers;
	//int f_num_init = params->num_f;
	//////////////////////////////////////

	int o_dim;
	int ch_num=1;//init
	struct conv_data_ *ptr_conv_data = &conv_data;
	//struct gn_data_ *ptr_gn_data = &gn_data;
	struct act_func_data_ *ptr_act_func_data = &act_func_data;
	struct maxpool_data_ *ptr_maxpool_data = &maxpool_data;
	struct concat_crop_data_ *ptr_concat_crop_data =&concat_crop_data;


	float ***conv_out;
	float ****conv_arr = (float ****)malloc(5*sizeof(float ***));//save each conv#_2 so we can concat later

	for(int curr_layer=1; curr_layer<6; curr_layer++)
	{
		///////////////////////////////      LAYER 1          /////////////////////////////////
		//1st convolution
		//(with zero padding ='same',so with stride =1 we get same dim as the input)
		//input: (1,dim,dim), filter :(16,1,3,3), output: (16,dim,dim)


		ptr_conv_data->ch_num = ch_num;
		ptr_conv_data->dim = dim;
		ptr_conv_data->f_dim =3;
		ptr_conv_data->f_num = calc_f_num(curr_layer); //that will be the output channels
		ptr_conv_data->mode = 1; // padding = 'same'
		ptr_conv_data->conv_in = image;
		ptr_conv_data->filter = filters[(curr_layer-1)*2];//2*9 + 1 filters
		ptr_conv_data->bias = bias[(curr_layer-1)*2];	   //2*9 + 1 bias

		conv(ptr_conv_data);

		conv_out = ptr_conv_data->conv_out;
		dim = ptr_conv_data->o_dim;
		ch_num = ptr_conv_data->f_num;

		/*
		//GN
		ptr_gn_data->batch = gn_batch;
		ptr_gn_data->ch_num = ch_num;
		ptr_gn_data->dim = dim;
		ptr_gn_data->gamma=gamma[(curr_layer-1)*2];//2*9
		ptr_gn_data->beta= beta[(curr_layer-1)*2]; //2*9
		ptr_gn_data->image = conv_out;

		GN(ptr_gn_data);

		conv_out = ptr_gn_data->out;
		*/

		//Relu
		ptr_act_func_data->Z = conv_out;
		ptr_act_func_data->channels=ch_num;
		ptr_act_func_data->dim=dim;
		ptr_act_func_data->code=3; //code for relu==3


		Activation_Function(ptr_act_func_data);
		conv_out = ptr_act_func_data->res;
		//2nd convolution
		//input: (16,dim,dim), filter (16,16,3,3), output: (16,dim,dim)
		//conv_block()
		ptr_conv_data->ch_num = ch_num;// or calc_ch_num(layer,#conv);
		ptr_conv_data->dim = dim;
		ptr_conv_data->f_dim =3;
		ptr_conv_data->f_num = calc_f_num(curr_layer); //that will be the output channels
		ptr_conv_data->mode = 1; // padding = 'same'
		ptr_conv_data->conv_in = conv_out;
		ptr_conv_data->filter = filters[(curr_layer-1)*2+1];//2*9 + 1 filters
		ptr_conv_data->bias = bias[(curr_layer-1)*2+1];	   //2*9 + 1 bias

		conv(ptr_conv_data);

		conv_out = ptr_conv_data->conv_out;
		dim = ptr_conv_data->o_dim;
		ch_num = ptr_conv_data->f_num;

		/*
		//GN
		ptr_gn_data->batch = gn_batch;
		ptr_gn_data->ch_num = ch_num;
		ptr_gn_data->dim = dim;
		ptr_gn_data->gamma=gamma[(curr_layer-1)*2+1];//2*9
		ptr_gn_data->beta= beta[(curr_layer-1)*2+1]; //2*9
		ptr_gn_data->image = conv_out;

		GN(ptr_gn_data);

		conv_out = ptr_gn_data->out;
		*/
		//Relu
		ptr_act_func_data->Z = conv_out;
		ptr_act_func_data->channels=ch_num;
		ptr_act_func_data->dim=dim;
		ptr_act_func_data->code=3; //code for relu==3



		Activation_Function(ptr_act_func_data);
		conv_out = ptr_act_func_data->res;

		///////////////////////////////////////////////////////////////////////////////////////
		conv_arr[curr_layer -1] = conv_out;
		//////////// maxpool /////////
		if(curr_layer != 5)
		{
			ptr_maxpool_data->channels=ch_num;
			ptr_maxpool_data->dim=dim;
			ptr_maxpool_data->image = conv_out;
			maxpool(ptr_maxpool_data);
			dim=ptr_maxpool_data->o_dim;
			image=ptr_maxpool_data->output;
		}
	}





	////////////////  Transposed convolution  //////////////////////////
	//conv_out is the input image

	for(int curr_layer=6; curr_layer<10; curr_layer++)
	{
		ptr_conv_data->ch_num = ch_num;
		ptr_conv_data->dim = dim;
		ptr_conv_data->f_dim=2;
		ptr_conv_data->f_num=calc_f_num(curr_layer);
		ptr_conv_data->mode = 0; //
		ptr_conv_data->conv_in = conv_out;
		ptr_conv_data->filter = f_dc[curr_layer-6];
		ptr_conv_data->bias = b_dc[curr_layer-6];

		convTransp(ptr_conv_data);

		dim = ptr_conv_data->o_dim;
		conv_out = ptr_conv_data->conv_out;
		ch_num = ptr_conv_data->f_num;


		//Now its time for concatination with respect to channels
		ptr_concat_crop_data->image1 = conv_arr[9-curr_layer];
		ptr_concat_crop_data->image2 = conv_out;
		ptr_concat_crop_data->dim=dim;
		ptr_concat_crop_data->ch_num=ch_num;

		concat(ptr_concat_crop_data);

		ch_num = ptr_concat_crop_data->o_ch_num;
		conv_out = ptr_concat_crop_data->image3;



		// 2 convolutions(without Group Normalization#TODO)

		ptr_conv_data->ch_num = ch_num;
		ptr_conv_data->dim = dim;
		ptr_conv_data->f_dim =3;
		ptr_conv_data->f_num = calc_f_num(curr_layer); //that will be the output channels
		ptr_conv_data->mode = 1; // padding = 'same'
		ptr_conv_data->conv_in = conv_out;
		ptr_conv_data->filter = filters[(curr_layer-1)*2];//2*9 + 1 filters
		ptr_conv_data->bias = bias[(curr_layer-1)*2];	   //2*9 + 1 bias

		conv(ptr_conv_data);

		conv_out = ptr_conv_data->conv_out;
		dim = ptr_conv_data->o_dim;
		ch_num = ptr_conv_data->f_num;

		/*
		//GN
		ptr_gn_data->batch = gn_batch;
		ptr_gn_data->ch_num = ch_num;
		ptr_gn_data->dim = dim;
		ptr_gn_data->gamma=gamma[(curr_layer-1)*2];//2*9
		ptr_gn_data->beta= beta[(curr_layer-1)*2]; //2*9
		ptr_gn_data->image = conv_out;

		GN(ptr_gn_data);

		conv_out = ptr_gn_data->out;
		*/

		//Relu
		ptr_act_func_data->Z = conv_out;
		ptr_act_func_data->channels=ch_num;
		ptr_act_func_data->dim=dim;
		ptr_act_func_data->code=3; //code for relu==3

		Activation_Function(ptr_act_func_data);
		conv_out = ptr_act_func_data->res;


		//2nd convolution
		ptr_conv_data->ch_num = ch_num;// or calc_ch_num(layer,#conv);
		ptr_conv_data->dim = dim;
		ptr_conv_data->f_dim =3;
		ptr_conv_data->f_num = calc_f_num(curr_layer); //that will be the output channels
		ptr_conv_data->mode = 1; // padding = 'same'
		ptr_conv_data->conv_in = conv_out;
		ptr_conv_data->filter = filters[(curr_layer-1)*2+1];//2*9 + 1 filters
		ptr_conv_data->bias = bias[(curr_layer-1)*2+1];	   //2*9 + 1 bias

		conv(ptr_conv_data);

		conv_out = ptr_conv_data->conv_out;
		dim = ptr_conv_data->o_dim;
		ch_num = ptr_conv_data->f_num;

		/*
		//GN
		ptr_gn_data->batch = gn_batch;
		ptr_gn_data->ch_num = ch_num;
		ptr_gn_data->dim = dim;
		ptr_gn_data->gamma=gamma[(curr_layer-1)*2+1];//2*9
		ptr_gn_data->beta= beta[(curr_layer-1)*2+1]; //2*9
		ptr_gn_data->image = conv_out;

		GN(ptr_gn_data);

		conv_out = ptr_gn_data->out;
		*/

		//Relu
		ptr_act_func_data->Z = conv_out;
		ptr_act_func_data->channels=ch_num;
		ptr_act_func_data->dim=dim;
		ptr_act_func_data->code=3; //code for relu==3

		Activation_Function(ptr_act_func_data);
		conv_out = ptr_act_func_data->res;

	}


	////////// Last(single conv) layer !!!!! ////
	int curr_layer = 10;
	ptr_conv_data->ch_num = ch_num;// or calc_ch_num(layer,#conv);
	ptr_conv_data->dim = dim;
	ptr_conv_data->f_dim =1;
	ptr_conv_data->f_num = calc_f_num(curr_layer); //that will be the output channels
	ptr_conv_data->mode = 0; //
	ptr_conv_data->conv_in = conv_out;
	ptr_conv_data->filter = filters[(curr_layer-1)*2];//2*9 + 1 filters
	ptr_conv_data->bias = bias[(curr_layer-1)*2];	   //2*9 + 1 bias

	conv(ptr_conv_data);

	conv_out = ptr_conv_data->conv_out;
	dim = ptr_conv_data->o_dim;
	ch_num = ptr_conv_data->f_num;

	///Normalization for possible inf edges
	struct norm_data_ *ptr_norm_data=&norm_data;
	//ptr_norm_data->image=conv_out;
	//ptr_norm_data->dim=dim;
	//ptr_norm_data->code = 0;

	//normalize_custom(ptr_norm_data);//it changes the conv_out ifself

	//Relu
	ptr_act_func_data->Z = conv_out;
	ptr_act_func_data->channels=ch_num;
	ptr_act_func_data->dim=dim;
	ptr_act_func_data->code=1; //code for sigmoid(==1)

	Activation_Function(ptr_act_func_data);
	conv_out = ptr_act_func_data->res;


	// fix some values that didnt reach their target value
	ptr_norm_data->image=conv_out;
	ptr_norm_data->dim=dim;
	ptr_norm_data->code = 1;

	normalize_custom(ptr_norm_data);//it changes the conv_out itself

	float accuracy = Dice_Coef(conv_out, label,dim);
	printf("\n\nAccuracy: %.2f % \n\n", (accuracy*100));
}
