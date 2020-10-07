#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <main.h>
#include <dirent.h>
#include <string.h>

void load_images(struct images_data_ *images_data)
{

	char **image_names;
	int im_num = images_data->im_num;//NUmber of images in directory
	int dim = images_data->dim;     //Resolution of the images
	//create char space for image names
	image_names = (char **)malloc(im_num*sizeof(char *));
	for (int i = 0; i<im_num ; i++)
		image_names[i]=(char *)malloc(50*sizeof(char));

	DIR *dir;
	struct dirent *ent;

	//User must edit this path !
	if ((dir = opendir ("/home/labis/eclipse-workspace/Utilities/images")) != NULL)
	{
		printf("\nLoading Images . . .");
	    /* print all the files and directories within directory */
		int i=0;
		int count=0;
		while ((ent = readdir (dir)) != NULL)
		{
			//printf ("\n%s , type: %d\n", ent->d_name, ent->d_type);
			if(ent->d_type == 8)//shows that its a string name, then save it
			{
				strcpy(image_names[i], ent->d_name);
				printf("\n%d) %s\n",count++,image_names[i]);
				i++;
			}
		}
	  closedir (dir);
	}
	else
	{
	  /* could not open directory */
	  perror ("Couldnt Open Directory(Check path!)");
	  return exit(1);
	}

    char line[20], path_name[200];
    int width, height, maxval;
    float ****image;
    FILE *fd;
    int ch_num=1;
    image = make_4darray(im_num, ch_num, dim);//create the required space for image storing
    for(int im=0; im< im_num; im++)
    {
    	// PATH EDIT REQUIRED!!!
    	sprintf(path_name,"/home/labis/eclipse-workspace/Utilities/images/%s", image_names[im]);
    	fd = fopen(path_name,"rb");
		if (fd == NULL) {
			printf("Could not open image file!\n");
			exit(1);
		}
	    fgets(line, sizeof(line), fd);
	    if (strcmp(line, "P5\n") != 0) { // Read the 1st line and check if the file format is valid(.PGM)
	        printf("Image is not in PGM(P5) format!\n");
	        fclose(fd);
	        exit(1);
	    }


	    // read header(2nd line),(includes the widht and height resolution, there is a space character between them)
	    if (fscanf(fd, "%d %d\n", &width, &height) != 2) {
	        printf("Invalid header(width & height area)\n");
	        fclose(fd);
	        exit(1);
	    }
	    if (fscanf(fd, "%d\n", &maxval) != 1) {//3rd line: maxvalue: if 16bit-->65535, if 8 bit->256
	        printf("Invalid header(maxval area)\n");
	        fclose(fd);
	        exit(1);
	    }
	    if (maxval == 65535) {//then its a 16bit image

	    	int ch_num=1,dim=height;
			uint16_t *rbuffer = (uint16_t *)malloc(dim*dim*sizeof(int16_t));//create the appropriate 16bit buffer for the whole image
			if (fread(rbuffer, sizeof(uint16_t), (width * height), fd) != width * height) {//read all the image at once
				printf("Error reading pixel values!\n");
				fclose(fd);
				exit(1);
			}
			fclose(fd);
			int offset=0;
			for(int i = 0; i<ch_num; i++ )
			{
				for (int x=0; x<height; x++)
				{
					for (int y=0; y<width; y++)
					{
						image[im][0][x][y] = (float)((*(rbuffer+offset))/maxval);//normalized
						offset++;
						//printf("%.3f\t", image[i][x][y]);
					}
					//printf("\n");
				}
				//printf("\n");
			}
			//free(rbuffer);
	    }
	    else //means its 8-bit
	    {
	    	int ch_num=1,dim=height;
			uint8_t *rbuffer = (uint8_t *)malloc(dim*dim*sizeof(int8_t));
			if (fread(rbuffer, sizeof(uint8_t), (width * height), fd) != width * height) {
				printf("Error reading pixel values!\n");
				fclose(fd);
				exit(1);
			}
			fclose(fd);
			int offset=0;
			for(int i = 0; i<ch_num; i++ )
			{
				for (int x=0; x<height; x++)
				{
					for (int y=0; y<width; y++)
					{
						image[im][0][x][y] = ((float)*(rbuffer+offset))/maxval;//normalized
						offset++;
						//printf("%.3f\t", image[i][x][y]);
					}
					//printf("\n");
				}
				//printf("\n");
			}
			//free(rbuffer);
	    }
    }
    images_data->images = image;
    printf("Done!\n");

}

//Same for labels
void load_labels(struct images_data_ *images_data)
{

	char **label_names;
	int im_num = images_data->im_num;
	int dim = images_data->dim;
	//create char space forr image names
	label_names = (char **)malloc(im_num*sizeof(char *));
	for (int i = 0; i<im_num ; i++)
		label_names[i]=(char *)malloc(50*sizeof(char));

	DIR *dir;
	struct dirent *ent;
	if ((dir = opendir ("/home/labis/eclipse-workspace/Utilities/labels")) != NULL) {
		printf("Loading Labels . . .");
	  /* print all the files and directories within directory */
		int i=0;
		while ((ent = readdir (dir)) != NULL)
		{
	    //printf ("%s\n", ent->d_name);
			if(ent->d_type == 8)//shows that its a string name
			{

				strcpy(label_names[i], ent->d_name);
				//printf("\n%s\n",label_names[i]);
				i++;
			}
		}
	  closedir (dir);
	}
	else
	{
	  /* could not open directory */
	  perror ("Coudnt open directory(Check path!)");
	  return exit(1);
	}
	//"/home/labis/data/salt/testfile.bin","rb"
    char line[20], path_name[200];
    int width, height, maxval;
    float ****label;
    FILE *fd;
    int ch_num=1;
    label = make_4darray(im_num, ch_num, dim);
    for(int im=0; im< im_num; im++)
    {
    	sprintf(path_name,"/home/labis/eclipse-workspace/Utilities/labels/%s", label_names[im]);
    	fd = fopen(path_name,"rb");
		if (fd == NULL) {
			printf("Could not open image file!\n");
			exit(1);
		}
	    fgets(line, sizeof(line), fd);
	    if (strcmp(line, "P5\n") != 0) {
	        printf("Image is not in PGM(P5) format!\n");
	        fclose(fd);
	        exit(1);
	    }

	    // skip comment
	    //fgets(line, sizeof(line), fd);

	    // read header
	    if (fscanf(fd, "%d %d\n", &width, &height) != 2) {
	        printf("Invalid header(width & height area)\n");
	        fclose(fd);
	        exit(1);
	    }
	    if (fscanf(fd, "%d\n", &maxval) != 1) {
	        printf("Invalid header(maxval area)\n");
	        fclose(fd);
	        exit(1);
	    }
	    if (maxval == 65535) {

	    	int ch_num=1,dim=height;
			uint16_t *rbuffer = (uint16_t *)malloc(dim*dim*sizeof(int16_t));
			if (fread(rbuffer, sizeof(uint16_t), (width * height), fd) != width * height) {
				printf("Error reading pixel values!\n");
				fclose(fd);
				exit(1);
			}
			fclose(fd);
			int offset=0;
			for(int i = 0; i<ch_num; i++ )
			{
				for (int x=0; x<height; x++)
				{
					for (int y=0; y<width; y++)
					{
						label[im][0][x][y] = (float)((*(rbuffer+offset))/maxval);//normalized
						offset++;
						//printf("%.3f\t", image[i][x][y]);
					}
					//printf("\n");
				}
				//printf("\n");
			}
	    }
	    else //means its 8-bit
	    {
	    	int ch_num=1,dim=height;
			uint8_t *rbuffer = (uint8_t *)malloc(dim*dim*sizeof(int8_t));
			if (fread(rbuffer, sizeof(uint8_t), (width * height), fd) != width * height) {
				printf("Error reading pixel values!\n");
				fclose(fd);
				exit(1);
			}
			fclose(fd);
			int offset=0;
			for(int i = 0; i<ch_num; i++ )
			{
				for (int x=0; x<height; x++)
				{
					for (int y=0; y<width; y++)
					{
						label[im][0][x][y] = ((float)*(rbuffer+offset))/maxval;//normalized
						offset++;
						//printf("%.3f\t", image[i][x][y]);
					}
					//printf("\n");
				}
				//printf("\n");
			}
	    }
    }
    images_data->labels = label;
    printf("Done!\n");
}

void load_params(struct params_ *params)
{
	// Unpacking //

	int batch =params->gn_batch; //      default:2
	int layers = params->layers; //  !!! default:10  !!!
	//int f_num = params->num_f;
	/////////////////
	float *****filters = (float *****)malloc((layers*2*2-1)*sizeof(float ****));
	float **bias = (float **)malloc((layers*2*2-1)*sizeof(float *));
	float *****f_dc = (float *****)malloc((layers-2)*sizeof(float ****));//NO OUT layer
	float **b_dc = (float **)malloc((layers-2)*sizeof(float *));
	float **gamma = (float **)malloc((layers*2*2-2)*sizeof(float *));//no out layers
	float **beta = (float **)malloc((layers*2*2-2)*sizeof(float *));


	////////////////////////// FILTERS ////////////////////////////
	///////////////////////////////////////////////////////////////
	int dim=3;
	int f_num,ch_num;
	int sum=0;
	int offset=0;

	for (int i=1;i<=layers; i++)// layers : 10 (default)
	{
		if(i!=10) // last layer is a special 1x1 convolution
		{
			sum += calc_f_num(i)*calc_ch_num(i,1)*3*3;//calculating number of parameters (1st part of the convolution block)
			sum += calc_f_num(i)*calc_ch_num(i,2)*3*3;//calculating number of parameters (2nd part of the convolution block)
		}
		else
		{
			sum+=calc_f_num(i)*calc_ch_num(i,1)*1*1; //last layer: out_f
		}
	}
	uint32_t *rbuffer;
	FILE *ptr = fopen("/home/labis/eclipse-workspace/Utilities/weights_encrypted.bin","rb"); // EDIT PATH //
	if(ptr == NULL)
	{
		printf("Couldnt load directory!\nExiting . . . ");
		exit(1);
	}
	printf("Loading Parameters . . .");
	rbuffer = (uint32_t *)malloc(sum*sizeof(int32_t));
	fread(rbuffer, sum*sizeof(uint32_t), 1, ptr);
	int pos=0;
	for (int i=1;i<=layers; i++)
	{
		if(i!=10)// All layer except last one which is the 1x1 convolution
		{
			f_num = calc_f_num(i);    //calculatre current layer's filters
			ch_num = calc_ch_num(i,1);//calculate current layer's channels(filter)
			float ****f = make_4darray(calc_f_num(i), calc_ch_num(i,1), 3);
			for(int k=0; k<f_num; k++)      //filters
				for(int j=0; j< ch_num; j++)//channels
					for(int x=0; x<dim; x++)
						for (int y=0; y<dim; y++)
						{
							f[k][j][x][y] = *((float *)(rbuffer+offset));
							offset++;
						}
			pos = (i-1)*2;
			filters[pos]=f; // careful - specific way of indexing data

			//f_num = calc_f_num(i);
			ch_num = calc_ch_num(i,2);
			f = make_4darray(f_num, ch_num, 3);
			for(int k=0; k<f_num; k++)
				for(int j=0; j< ch_num; j++)
					for(int x=0; x<dim; x++)
						for (int y=0; y<dim; y++)
						{
							f[k][j][x][y] = *((float *)(rbuffer+offset));
							offset++;
						}
			pos += 1;
			filters[pos]=f;
		}
		else//i == 10 --> last 1x1 convolution
		{
			f_num = calc_f_num(i);
			ch_num = calc_ch_num(i,1);
			float ****f = make_4darray(f_num, ch_num, 1);
			for(int k=0; k<f_num; k++)
				for(int j=0; j< ch_num; j++)
					for(int x=0; x<1; x++)
						for (int y=0; y<1; y++)
						{
							f[k][j][x][y] = *((float *)(rbuffer+offset));
							offset++;
						}
			pos = (i-1)*2; //pos == 18 (19 filters)
			filters[pos]=f;
		}
	}

	/////////////////// BIAS ///////////////////
	////////////////////////////////////////////
	sum=1;//already last layer out is counted, (just 1 scalar bias_out)
	offset=0;
	for (int i=1;i<=(layers-1); i++)//9 layers, last layer is calculated above(sum=1)
		sum += calc_f_num(i)*2;

	rbuffer = (uint32_t *)malloc(sum*sizeof(int32_t));
	fread(rbuffer, sum*sizeof(uint32_t), 1, ptr);
	pos=0;
	offset=0;
	for (int i=1;i<=(layers-1); i++) //NOT THE LAST LAYER!!!
	{
		f_num = calc_f_num(i);
		ch_num = calc_ch_num(i,1);
		float *b = (float *)malloc(f_num*sizeof(float));
		for (int y=0; y<f_num; y++)
		{
			b[y] = *((float *)(rbuffer+offset));
			offset++;
		}
		pos = (i-1)*2;
		bias[pos]=b;

		//f_num = calc_f_num(i);
		ch_num = calc_ch_num(i,2);
		b =  (float *)malloc(f_num*sizeof(float));
		for (int y=0; y<f_num; y++)
		{
			b[y] = *((float *)(rbuffer+offset));
			offset++;
		}
		pos += 1;
		bias[pos]=b;
	}
	//last layer of bias
	f_num = calc_f_num(10);
	float *b = (float *)malloc(f_num*sizeof(float));
	for (int y=0; y<f_num; y++)
	{
		b[y] = *((float *)(rbuffer+offset));
		offset++;
	}
	bias[18]=b;

	///////////////////////////// F_DC ///////////////////////////////
	//////////////////////////////////////////////////////////////////

	dim=2;
	sum=0;
	offset=0;
	sum =((128*256) +(64*128)+(32*64)+(16*32))*2*2;// calculation: filter1_num*channel1_num + ... + filterN_numb*channelN_num)*KernelSize*KernelSize
	rbuffer = (uint32_t *)malloc(sum*sizeof(int32_t));
	fread(rbuffer, sum*sizeof(uint32_t), 1, ptr);
	pos=0;
	for (int i=6;i<=(layers-1); i++)
	{
		f_num = calc_f_num(i);//same as filter i_1
		ch_num = calc_ch_num(i,1);//same as filter i_1 (e.g. Transp conv 6 filter =(shape)= conv 6_1)
		float ****f = make_4darray(calc_f_num(i), calc_ch_num(i,1), 2);
		for(int k=0; k<f_num; k++)
			for(int j=0; j< ch_num; j++)
				for(int x=0; x<dim; x++)
					for (int y=0; y<dim; y++)
					{
						f[k][j][x][y] = *((float *)(rbuffer+offset));
						offset++;
					}
		f_dc[pos]=f;
		pos++;
	}
	///////////////////////////// B_DC ///////////////////////////////
	//////////////////////////////////////////////////////////////////

	sum=0;
	offset=0;
	sum = 128+64+32+16;
	rbuffer = (uint32_t *)malloc(sum*sizeof(int32_t));
	fread(rbuffer, sum*sizeof(uint32_t), 1, ptr);
	pos=0;
	for (int i=6;i<=(layers-1); i++)
	{
		f_num = calc_f_num(i);//same as filter 6_1
		float *b = (float *)malloc(f_num*sizeof(float));
		for (int y=0; y<f_num; y++)
		{
			b[y] = *((float *)(rbuffer+offset));
			offset++;
		}
		b_dc[pos]=b;
		pos++;
	}
	/*
	/////////////////// GAMMA ///////////////////
	////////////////////////////////////////////
	sum=0;//already last layer out is counted, (just 1 scalar bias_out)
	offset=0;
	for (int i=1;i<=(layers-1); i++)//9 layers
		sum += ((int)(calc_f_num(i)/batch))*2;//default:batch=2

	rbuffer = (uint32_t *)malloc(sum*sizeof(int32_t));
	fread(rbuffer, sum*sizeof(uint32_t), 1, ptr);
	pos=0;
	offset=0;
	for (int i=1;i<=(layers-1); i++) //NOT THE LAST LAYER!!!
	{
		f_num = (int)(calc_f_num(i)/batch);
		float *ga = (float *)malloc(f_num*sizeof(float));
		for (int y=0; y<f_num; y++)
		{
			ga[y] = *((float *)(rbuffer+offset));
			offset++;
		}
		pos = (i-1)*2;
		gamma[pos]=ga;


		ga =  (float *)malloc(f_num*sizeof(float));
		for (int y=0; y<f_num; y++)
		{
			ga[y] = *((float *)(rbuffer+offset));
			offset++;
		}
		pos += 1;
		gamma[pos]=ga;
	}
	/////////////////// BETA ///////////////////
	////////////////////////////////////////////
	sum=0;//already last layer out is counted, (just 1 scalar bias_out)
	offset=0;
	for (int i=1;i<=(layers-1); i++)//9 layers
		sum += ((int)(calc_f_num(i)/batch))*2;//default:batch=2

	rbuffer = (uint32_t *)malloc(sum*sizeof(int32_t));
	fread(rbuffer, sum*sizeof(uint32_t), 1, ptr);
	pos=0;
	offset=0;
	for (int i=1;i<=(layers-1); i++) //NOT THE LAST LAYER!!!
	{
		f_num = (int)(calc_f_num(i)/batch);
		float *be = (float *)malloc(f_num*sizeof(float));
		for (int y=0; y<f_num; y++)
		{
			be[y] = *((float *)(rbuffer+offset));
			offset++;
		}
		pos = (i-1)*2;
		beta[pos]=be;


		be =  (float *)malloc(f_num*sizeof(float));
		for (int y=0; y<f_num; y++)
		{
			be[y] = *((float *)(rbuffer+offset));
			offset++;
		}
		pos += 1;
		beta[pos]=be;
	}
	params->be=beta;
	params->ga=gamma;
	*/

	//Pack parameters back to the structure
	params->b_dc=b_dc;
	params->f_dc = f_dc;
	params->bias=bias;
	params->filters=filters;
	fclose(ptr);
	printf("Done!\n");

}
