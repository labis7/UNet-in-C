#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <main.h>


void load_params(struct params_ *params)
{
	// Unpacking //

	int batch =params->gn_batch; //default:2
	int layers = params->layers; //default:10
	//int f_num = params->num_f;
	/////////////////
	float *****filters = (float *****)malloc((layers*2*2-1)*sizeof(float ****));
	float **bias = (float **)malloc((layers*2*2-1)*sizeof(float *));
	float *****f_dc = (float *****)malloc((layers-2)*sizeof(float ****));//NO OUT layer
	float **b_dc = (float **)malloc((layers-2)*sizeof(float *));
	float **gamma = (float **)malloc((layers*2*2-2)*sizeof(float *));//no out layers
	float **beta = (float **)malloc((layers*2*2-2)*sizeof(float *));

	//init reading from file //

	////////////////////////// FILTERS ////////////////////////////
	///////////////////////////////////////////////////////////////
	int dim=3;
	int f_num,ch_num;
	int sum=0;
	int offset=0;
	for (int i=1;i<=layers; i++)
	{
		if(i!=10)
		{
			sum += calc_f_num(i)*calc_ch_num(i,1)*3*3;
			sum += calc_f_num(i)*calc_ch_num(i,2)*3*3;
		}
		else
		{
			sum+=calc_f_num(i)*calc_ch_num(i,1)*1*1; //last layer: out_f
		}
	}
	uint32_t *rbuffer;
	FILE *ptr = fopen("/home/labis/data/salt/testfile.bin","rb");
	rbuffer = (uint32_t *)malloc(sum*sizeof(int32_t));
	fread(rbuffer, sum*sizeof(uint32_t), 1, ptr);
	int pos=0;
	for (int i=1;i<=layers; i++)
	{
		if(i!=10)
		{
			f_num = calc_f_num(i);
			ch_num = calc_ch_num(i,1);
			float ****f = make_4darray(calc_f_num(i), calc_ch_num(i,1), 3);
			for(int k=0; k<f_num; k++)
				for(int j=0; j< ch_num; j++)
					for(int x=0; x<dim; x++)
						for (int y=0; y<dim; y++)
						{
							f[k][j][x][y] = *((float *)(rbuffer+offset));
							offset++;
						}
			pos = (i-1)*2;
			filters[pos]=f;

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
		else//i == 10 --> last 1x1 conv
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
	for (int i=1;i<=(layers-1); i++)//9 layers, last layer is precalculated
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
	sum =((128*256) +(64*128)+(32*64)+(16*32))*2*2;
	rbuffer = (uint32_t *)malloc(sum*sizeof(int32_t));
	fread(rbuffer, sum*sizeof(uint32_t), 1, ptr);
	pos=0;
	for (int i=6;i<=(layers-1); i++)
	{
		f_num = calc_f_num(i);//same as filter 6_1
		ch_num = calc_ch_num(i,1);//same as filter 6_1
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
	params->b_dc=b_dc;
	params->f_dc = f_dc;
	params->bias=bias;
	params->filters=filters;
	fclose(ptr);

}

void GN(struct gn_data_ *gn_data)
{
	int eps=1e-5;
	float *beta,*gamma;
	float ivar, sqrtvar, var, mu;
	float ***xmu, ***xhat, ***gammax, ***out;
	int batch, ch_num, dim;
	batch=gn_data->batch;
	ch_num=gn_data->ch_num;
	dim=gn_data->dim;
	gamma=gn_data->gamma;
	beta=gn_data->beta;
	float ***image = gn_data->image;
	// make empty arrays that will need later
	xmu = make_3darray(batch,dim);//we need only size of batch in order to calculate the partial out on each loop
	xhat= make_3darray(batch,dim);
	gammax= make_3darray(batch,dim);
	out = make_3darray(ch_num,dim);
	/////////////////////////////////////////
	for (int i=0; i<ch_num; i+=batch)
	{
		//each loop is interested on i->i+batch channels,
		//with the data(gamma,beta) be the elements (int)(i//batch)

		//step1:calculate mean(scalar)
		mu=0;
		for(int k=i; k<(i+batch); k++)
			for(int x=0; x<dim; x++)
				for(int y=0; y<dim; y++)
					mu += image[k][x][y];
		mu = (float)(mu/(dim*dim*batch));

		//step2:subtract mean vector of every training example
		for(int k=i; k<(i+batch); k++)
			for(int x=0; x<dim; x++)
				for(int y=0; y<dim; y++)
					xmu[k-i][x][y] = image[k][x][y] - mu;

		//step3: following the lower branch - calculation denominator
		//step4: calculate variance
		var=0;
		for(int k=0; k<batch; k++)
			for(int x=0; x<dim; x++)
				for(int y=0; y<dim; y++)
					var+=pow(xmu[k][x][y],2);
		var = (float)(var/(batch*dim*dim));

		//step5: add eps for numerical stability, then sqrt
		sqrtvar = sqrt((var+eps));
		//step6: invert sqrtvar
		ivar = (float)(1/sqrtvar);

		//step7: execute normalization
		for(int k=0; k<batch; k++)
			for(int x=0; x<dim; x++)
				for(int y=0; y<dim; y++)
					xhat[k][x][y] = xmu[k][x][y] * ivar;


		//step8: Nor the two transformation steps
		int pos =(int)(i/batch);
		float gamma_t=gamma[pos];
		for(int k=0; k<batch; k++)
			for(int x=0; x<dim; x++)
				for(int y=0; y<dim; y++)
					gammax[k][x][y]=gamma_t*xhat[k][x][y];
		//step9: output
		float beta_t = beta[(int)(i/batch)];
		for(int k=i; k<(i+batch); k++)
			for(int x=0; x<dim; x++)
				for(int y=0; y<dim; y++)
					out[k][x][y]= gammax[k-i][x][y] + beta_t;
	}
	gn_data->out = out;
	//step5: add eps for numerical stability, then sqrt

}