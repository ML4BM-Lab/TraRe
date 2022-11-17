//variational Bayes spike regression (vbsr) R package C library
//Copyright 2012 Benjamin A Logsdon
#include "vbsr.h"


//xc: return jth column of data matrix x
double * xc(struct model_struct * model, int j){
	return (&(model->data.X[j]))->col;
}

double * xcm(struct model_marg_struct * model, int j){
	return (&(model->data.X[j]))->col;
}

//oc: return jth column of ordering matrix "ordering"
int * oc(struct model_struct * model, int j){
	return (&(model->data.ordering[j]))->col;
}

//me: return the ith, jth element of the ordering, path model parameters
struct model_param_struct * me(struct model_struct * model,int i, int j){
	return (&((&(model->order[i]))->model_param[j]));
}



void initialize_model_param(int n, 
				int m,
				int i, 
				int j,
				struct model_struct * model, 
				double * y, 
				double var_y){


	int k;

	me(model,i,j)->beta_mu = (double *) malloc(sizeof(double)*m);
	me(model,i,j)->beta_sigma = (double *) malloc(sizeof(double)*m);
	me(model,i,j)->beta_chi = (double *) malloc(sizeof(double)*m);
	me(model,i,j)->beta_p = (double *) malloc(sizeof(double)*m);
	me(model,i,j)->e_beta = (double *) malloc(sizeof(double)*m);
	me(model,i,j)->e_beta_sq = (double *) malloc(sizeof(double)*m);


	for(k=0;k<m;k++){
		me(model,i,j)->beta_mu[k] = 0;
		me(model,i,j)->beta_sigma[k] = 0;
		me(model,i,j)->beta_chi[k]= 0;
		me(model,i,j)->beta_p[k] = 0;
		me(model,i,j)->e_beta[k] = 0;
		me(model,i,j)->e_beta_sq[k] = 0;
	}


	me(model,i,j)->sigma_e = var_y;
	me(model,i,j)->lb = 0;
	me(model,i,j)->p_sums = 0;
	me(model,i,j)->entropy = 0;
	me(model,i,j)->v_sums_correct = 0;
	me(model,i,j)->w_vec = (double *) malloc(sizeof(double)*n);
	me(model,i,j)->mu_vec = (double *) malloc(sizeof(double)*n);
	me(model,i,j)->resid_vec = (double *) malloc(sizeof(double)*n);
	me(model,i,j)->pred_vec_old = (double *) malloc(sizeof(double)*n);
	me(model,i,j)->pred_vec_new = (double *) malloc(sizeof(double)*n);
	me(model,i,j)->x_w = (double *) malloc(sizeof(double)*n);

	for(k=0;k<n;k++){
		me(model,i,j)->w_vec[k] = 0.25;
		me(model,i,j)->mu_vec[k] = 0.5;
		switch(model->control_param.regressType){
			case LINEAR:
				me(model,i,j)->resid_vec[k] = y[k];
				break;
			case LOGISTIC:
				me(model,i,j)->resid_vec[k] = (y[k] - 0.5)/(0.25);
				break;
		}
		me(model,i,j)->pred_vec_old[k] = 0;
		me(model,i,j)->pred_vec_new[k] = 0;
		me(model,i,j)->x_w[k]=0;
	}

	me(model,i,j)->ord_index = i;
	me(model,i,j)->path_index = j;
}

void free_model_param(struct model_struct * model, int i, int j){

	free(me(model,i,j)->beta_mu);
	free(me(model,i,j)->beta_sigma);
	free(me(model,i,j)->beta_chi);
	free(me(model,i,j)->beta_p);
	free(me(model,i,j)->e_beta);
	free(me(model,i,j)->e_beta_sq);

	free(me(model,i,j)->w_vec);
	free(me(model,i,j)->mu_vec);
	free(me(model,i,j)->resid_vec);
	free(me(model,i,j)->pred_vec_old);
	free(me(model,i,j)->pred_vec_new);
	free(me(model,i,j)->x_w);

}


void ddot_w(int n,double *vect1,double *vect2,double * result){
	const int incxy = 1;
	(*result)=F77_NAME(ddot)(&n,vect1,&incxy,vect2,&incxy);
}


void daxpy_w(int n,double *x,double *y,double alpha){
	//y<- ax+y;
	const int incxy =1;
	F77_NAME(daxpy)(&n,&alpha,x,&incxy,y,&incxy);
}

void dnrm2_w(int n,double *x,double *result){
	const int incxy=1;
	(*result)=F77_NAME(dnrm2)(&n,x,&incxy);
}

void dscal_w(int n,double *x, double alpha){
	const int incxy=1;
	F77_NAME(dscal)(&n,&alpha,x,&incxy);
}


void scale_vector(double * vec,double * ones,int n){
	//mean zero, variance 1
	double mean,sd;
	double nd = (double) n;
	ddot_w(n,vec,ones,&mean);
	mean = mean/nd;
	daxpy_w(n,ones,vec,-mean);
	dnrm2_w(n,vec,&sd);
	dscal_w(n,vec,sqrt(nd)/(sd));
}


void cor(double * vec1, double * vec2, double * ones,double * corv,int n){

	double mean1,mean2,sd1,sd2,cov;
	double nd = (double) n;

	//compute means:	
	ddot_w(n,vec1,ones,&mean1);
	mean1 = mean1/nd;

	ddot_w(n,vec2,ones,&mean2);
	mean2 = mean2/nd;

	//compute variances:
	//substract means from vectors
	daxpy_w(n,ones,vec1,-mean1);
	dnrm2_w(n,vec1,&sd1);
	
	daxpy_w(n,ones,vec2,-mean2);
	dnrm2_w(n,vec2,&sd2);

	//rescale
	//sd1 = sd1/sq;
	//sd2 = sd2/sqrt(nd-1);

	ddot_w(n,vec1,vec2,&cov);
	//rescale
	//cov = cov/(;
	//Rprintf("cov: %g\n",cov);
	corv[0] = cov/(sd1*sd2);			

	//add means back to vectors
	daxpy_w(n,ones,vec1,mean1);
	daxpy_w(n,ones,vec2,mean2);


}

double compute_ssq(double *vec,int n){
	double a;
	ddot_w(n,vec,vec,&a);
	return a;
}


void process_data(struct model_struct * model){
	int j;
	double nd = ((double) model->data.n);

	switch(model->control_param.scaleType){

		case SCALE:
			//Rprintf("Scaling...\n");
			
			for(j=0;j<model->data.m;j++){
				if(j>0){
					scale_vector(xc(model,j),model->data.one_vec,model->data.n);
				}
				model->data.x_sum_sq[j] = nd;
			}
			break;

		case NOSCALE:
			//Rprintf("Sum of squares pre-compute...\n");
			for(j=0;j<model->data.m;j++){
				model->data.x_sum_sq[j]=compute_ssq(xc(model,j),model->data.n);
			}
			break;
	}


}


void process_data_marg(struct model_marg_struct * model){
	int j;
	double nd = ((double) model->data.n);

	switch(model->control_param.scaleType){

		case SCALE:
			//Rprintf("Scaling...\n");
			
			for(j=0;j<model->data.m;j++){
				if(j>0){
					scale_vector(xcm(model,j),model->data.one_vec,model->data.n);
				}
				model->data.x_sum_sq[j] = nd;
			}
			break;

		case NOSCALE:
			//Rprintf("Sum of squares pre-compute...\n");
			for(j=0;j<model->data.m;j++){
				model->data.x_sum_sq[j]=compute_ssq(xcm(model,j),model->data.n);
			}
			break;
	}


}


void initialize_model(double * eps, 
					double * l0_path, 
					double * pb_path, 
					int * exclude, 
					double * penalty_factor, 
					int * maxit, 
					int * path_length, 
					int * n_orderings, 
					int * regress, 
					int * scale, 
					int * est, 
					int * error, 
					double * kl,
					int * approx,
					int * total_replicates,
					double * X, 
					double * y, 
					double * var_y, 
					int * n, 
					int * m, 
					int * ordering_mat, 
					struct model_struct * model){

	//initialize: (*model).control_param;
	
	//Rprintf("eps: %g\n",eps[0]);
	//Rprintf("l0_path[0]: %g\n",l0_path[0]);
	//Rprintf("pb_path[0]: %g\n",pb_path[0]);
	//Rprintf("exclude[1]: %d\n",exclude[1]);
	//Rprintf("pf[0]: %g\n",penalty_factor[0]);
	//Rprintf("maxit: %d\n",maxit[0]);
	//Rprintf("path_length: %d\n",path_length[0]);
	//Rprintf("n_orderings: %d\n",n_orderings[0]);
	//Rprintf("regress: %d\n",regress[0]);
	//Rprintf("scale: %d\n",scale[0]);
	//Rprintf("est: %d\n",est[0]);
	//Rprintf("error: %d\n",error[0]);
	//Rprintf("kl: %g\n",kl[0]);
	//Rprintf("tr: %d\n",total_replicates[0]);
	//Rprintf("X[0,0]: %g\n",X[0]);
	//Rprintf("y[0]: %g\n",y[0]);
	//Rprintf("var_y: %g\n",var_y[0]);
	//Rprintf("n: %d\n",n[0]);
	//Rprintf("m: %d\n",m[0]);
	//Rprintf("ordering_mat[0]: %d\n",ordering_mat[0]);
	
	int k,l;
	model->control_param.eps = (*eps);
	//model->control_param.max_pb = (*max_pb);
	model->control_param.l0_path = (double *) malloc(sizeof(double)*(*path_length));
	model->control_param.pb_path = (double *) malloc(sizeof(double)*(*path_length));
	for(k=0;k<*path_length;k++){
		model->control_param.l0_path[k]=l0_path[k];
		model->control_param.pb_path[k]=pb_path[k];
	}
	
	model->control_param.exclude = (int *) malloc(sizeof(int)*(*m));
	model->control_param.penalty_factor = (double *) malloc(sizeof(double)*(*m));
	for(k=0;k<*m;k++){
		model->control_param.exclude[k]=exclude[k];
		model->control_param.penalty_factor[k]=penalty_factor[k];
	}
	
	model->control_param.maxit = (*maxit);
	model->control_param.path_length = (*path_length);
	model->control_param.n_orderings = (*n_orderings);
	if((*regress)==1){
		model->control_param.regressType = LINEAR;
	} else{
		model->control_param.regressType = LOGISTIC;
	}

	if((*scale)==1){
		model->control_param.scaleType = SCALE;
	} else{
		model->control_param.scaleType = NOSCALE;
	}

	if((*est)==1){
		model->control_param.estType = BMA;
	} else{
		model->control_param.estType = MAXIMAL;
	}

	if((*error)==1){
		model->control_param.errType = KL;
	} else{
		model->control_param.errType = NOKL;
	}
	
	if((*approx)==1){
		model->control_param.bType = APPR;
	} else{
		model->control_param.errType = EXACT;
	}


	model->control_param.kl_percentile = (*kl);
	model->control_param.total_replicates = (*total_replicates);
	//initialize: (*model).(data);
	//struct single_mod *single_mods;
	//single_mods= (single_mod *) malloc(sizeof(single_mod)*(n_order+1));

	model->data.X = (struct matrix_v *) malloc(sizeof(struct matrix_v)*(*m));
	for(k=0;k<(*m);k++){
		(&(model->data.X[k]))->col = (double *) malloc(sizeof(double)*(*n));
	}

	for(k=0;k<(*m);k++){
		for(l=0;l<(*n);l++){
			(&(model->data.X[k]))->col[l] = X[k*(*n)+l];
		}
	}
	
	model->data.y = y;
	model->data.var_y = (*var_y);
	model->data.n = (*n);
	model->data.m = (*m);
	int (pii) =0;
	for(k=0;k<(*m);k++){
		//Rprintf("exclude[%d]:%d\n",k,exclude[k]);
		if(exclude[k]==1){
			//Rprintf("worked\n");
			++(pii);
		}
	}
	model->data.p = (pii);
	//Rprintf("model->data.p = %d\n",model->data.p);
	model->data.x_sum_sq = (double *) malloc(sizeof(double)*(*m));


	model->data.ordering = (struct matrix_i *) malloc(sizeof(struct matrix_i)*(*n_orderings));
	for(k=0;k<(*n_orderings);k++){
		(&(model->data.ordering[k]))->col = (int *) malloc(sizeof(int)*(*m));
	}

	for(k=0;k<(*n_orderings);k++){
		for(l=0;l<(*m);l++){
			(&(model->data.ordering[k]))->col[l] = ordering_mat[k*(*m)+l];
		}
	}

	model->data.one_vec = (double *) malloc(sizeof(double)*(*n));
	for(k=0;k<(*n);k++){
		model->data.one_vec[k]= 1.0;
	}


	process_data(model);

	//initialize: (*model).me(model,i,j);

	model->order = (struct order_struct *) malloc(sizeof(struct order_struct)*(*n_orderings));
	for(k=0;k<(*n_orderings);k++){
		(&(model->order[k]))->model_param = (struct model_param_struct *) malloc(sizeof(struct model_param_struct)*(*path_length));
	}

	for(k=0;k<(*n_orderings);k++){
		for(l=0;l<(*path_length);l++){
			initialize_model_param((*n),(*m),k,l,model,y,*var_y);
		}
	}


}



void initialize_model_marg(double * eps,
			int * exclude,
			int * maxit,
			int * regress,
			int * scale,
			double * X,
			double * y,
			double * var_y,
			int * n,
			int * m,
			struct model_marg_struct * model){

	int k,l;
	
	//Rprintf("eps: %g\n",eps[0]);
	//Rprintf("exclude[1]: %d\n",exclude[1]);
	//Rprintf("maxit: %d\n",maxit[0]);
	//Rprintf("regress: %d\n",regress[0]);
	//Rprintf("scale: %d\n",scale[0]);
	//Rprintf("X[0,0]: %g\n",X[0]);
	//Rprintf("y[0]: %g\n",y[0]);
	//Rprintf("var_y: %g\n",var_y[0]);
	//Rprintf("n: %d\n",n[0]);
	//Rprintf("m: %d\n",m[0]);
	
			
	model->control_param.eps = (*eps);
	model->control_param.exclude = (int *) malloc(sizeof(int)*(*m));
	for(k=0;k<*m;k++){
		model->control_param.exclude[k]=exclude[k];
	}
	model->control_param.maxit = (*maxit);
	if((*regress)==1){
		model->control_param.regressType = LINEAR;
	} else{
		model->control_param.regressType = LOGISTIC;
	}

	if((*scale)==1){
		model->control_param.scaleType = SCALE;
	} else{
		model->control_param.scaleType = NOSCALE;
	}
	
	model->model_param.beta_mu = (double *) malloc(sizeof(double)*(*m));
	model->model_param.beta_sigma = (double *) malloc(sizeof(double)*(*m));
	model->model_param.beta_chi = (double *) malloc(sizeof(double)*(*m));
	model->model_param.beta_p = (double *) malloc(sizeof(double)*(*m));
	model->model_param.sigma_e = *var_y;
	
	
	for(k=0;k<(*m);k++){
		model->model_param.beta_mu[k] = 0;
		model->model_param.beta_sigma[k] = 0;
		model->model_param.beta_chi[k] = 0;
		model->model_param.beta_p[k] = 0;
	}

	model->model_param.w_vec = (double *) malloc(sizeof(double)*(*n));
	model->model_param.mu_vec = (double *) malloc(sizeof(double)*(*n));
	model->model_param.resid_vec = (double *) malloc(sizeof(double)*(*n));
	model->model_param.pred_vec_old = (double *) malloc(sizeof(double)*(*n));
	model->model_param.pred_vec_new = (double *) malloc(sizeof(double)*(*n));
	model->model_param.x_w = (double *) malloc(sizeof(double)*(*n));
  
	for(k=0;k<(*n);k++){
		model->model_param.w_vec[k] = 0.25;
		model->model_param.mu_vec[k] = 0.5;
		switch(model->control_param.regressType){
			case LINEAR:
				model->model_param.resid_vec[k] = y[k];
				break;
			case LOGISTIC:
				model->model_param.resid_vec[k] = (y[k] - 0.5)/(0.25);
				break;
		}
		model->model_param.pred_vec_old[k] = 0;
		model->model_param.pred_vec_new[k] = 0;
		model->model_param.x_w[k]=0;
	}
	
	
			
	model->data.X = (struct matrix_v *) malloc(sizeof(struct matrix_v)*(*m));
	for(k=0;k<(*m);k++){
		(&(model->data.X[k]))->col = (double *) malloc(sizeof(double)*(*n));
	}

	for(k=0;k<(*m);k++){
		for(l=0;l<(*n);l++){
			(&(model->data.X[k]))->col[l] = X[k*(*n)+l];
		}
	}
	
	model->data.y = y;
	model->data.var_y = (*var_y);
	model->data.n = (*n);
	model->data.m = (*m);
	model->data.x_sum_sq = (double *) malloc(sizeof(double)*(*m));
			
	model->data.one_vec = (double *) malloc(sizeof(double)*(*n));
	for(k=0;k<(*n);k++){
		model->data.one_vec[k]= 1.0;
	}


	process_data_marg(model);			
			
			
}





void free_model(struct model_struct * model){
	//free X
	int i,j,k;
	for(k=0;k<(model->data.m);k++){
		free((&(model->data.X[k]))->col);
		
	}
	free(model->data.X);
	//free orderings

	for(k=0;k<(model->control_param.n_orderings);k++){
		free((&(model->data.ordering[k]))->col);
	}
	free(model->data.ordering);
	
	//free order:model_param
	for(i=0;i<model->control_param.n_orderings;i++){
		for(j=0;j<model->control_param.path_length;j++){
			free_model_param(model, i, j);
		}
	}
	

	for(k=0;k<(model->control_param.n_orderings);k++){
		free((&(model->order[k]))->model_param);
	}
	free(model->order);
	
	//free x_sum_sq
	
	free(model->data.x_sum_sq);

	//free one_vec
	
	free(model->data.one_vec);
	
	free(model->control_param.l0_path);
	free(model->control_param.pb_path);
	free(model->control_param.penalty_factor);
	free(model->control_param.exclude);


}

void free_model_marg(struct model_marg_struct * model){

	int k;

	for(k=0;k<(model->data.m);k++){
		free((&(model->data.X[k]))->col);
		
	}
	free(model->data.X);

	free(model->model_param.beta_mu);
	free(model->model_param.beta_sigma);
	free(model->model_param.beta_chi);
	free(model->model_param.beta_p);

	free(model->model_param.w_vec);
	free(model->model_param.mu_vec);
	free(model->model_param.resid_vec);
	free(model->model_param.pred_vec_old);
	free(model->model_param.pred_vec_new);
	free(model->model_param.x_w);
	
	free(model->data.x_sum_sq);
	free(model->data.one_vec);
	free(model->control_param.exclude);
	

}



void copy_model_state(struct model_struct * model, int i, int j){
	int k,l;
	l = j-1;

	for(k=0;k<model->data.m;k++){
		me(model,i,j)->beta_mu[k] = me(model,i,l)->beta_mu[k];
		me(model,i,j)->beta_sigma[k] = me(model,i,l)->beta_sigma[k];
		me(model,i,j)->beta_chi[k] = me(model,i,l)->beta_chi[k];
		me(model,i,j)->beta_p[k] = me(model,i,l)->beta_p[k];
		me(model,i,j)->e_beta[k] = me(model,i,l)->e_beta[k];
		me(model,i,j)->e_beta_sq[k] = me(model,i,l)->e_beta_sq[k];
	}


	me(model,i,j)->sigma_e = me(model,i,l)->sigma_e;
	me(model,i,j)->lb = me(model,i,l)->lb;
	me(model,i,j)->p_sums = me(model,i,l)->p_sums;
	me(model,i,j)->entropy = me(model,i,l)->entropy;
	me(model,i,j)->v_sums_correct = me(model,i,l)->v_sums_correct;

	for(k=0;k<model->data.n;k++){
		me(model,i,j)->w_vec[k] = me(model,i,l)->w_vec[k];
		me(model,i,j)->mu_vec[k] = me(model,i,l)->mu_vec[k];
		me(model,i,j)->resid_vec[k] = me(model,i,l)->resid_vec[k];
		me(model,i,j)->pred_vec_old[k] = me(model,i,l)->pred_vec_old[k];
		me(model,i,j)->pred_vec_new[k] = me(model,i,l)->pred_vec_new[k];
	}



}



void update_beta(struct model_struct * model, int i, int j){

	int k,l,exc,t;
	double mu, sigma,prec, chi, p, e_b,e_b2,l0;


	//if(model->control_param.max_pb==1){
	//	l0 = me(model,i,j)->l0_max;
	//}else{
		l0 = model->control_param.l0_path[j];
	//}
	//Rprintf("l0: %g\n",l0);
	//Rprintf("m: %d\n",model->data.m);
	//error("!!!\n");
	switch (model->control_param.regressType){


		case LINEAR:
			//run linear updates
			for(l=0;l< model->data.m ;l++){
				k = (&(model->data.ordering[i]))->col[l];
				//k = l;
				//Rprintf("k: %d\n",k);
				exc = model->control_param.exclude[k];
				//Rprintf("exc: %d\n",exc);
				ddot_w(model->data.n,xc(model,k),me(model,i,j)->resid_vec,&mu);
				//Rprintf("mu: %g\n",mu);
				mu = mu + (model->data.x_sum_sq[k])*(me(model,i,j)->e_beta[k]);
				//Rprintf("mu: %g\n",mu);
				mu = mu/model->data.x_sum_sq[k];
				//Rprintf("mu: %g\n",mu);
				sigma = 1/((1/me(model,i,j)->sigma_e)*(model->data.x_sum_sq[k]));
				//Rprintf("sigma: %g\n",sigma);
				chi = pow(mu,2)/sigma;
				//Rprintf("chi: %g, e_beta[%d]: %g\n",chi,k,me(model,i,j)->e_beta[k]);
				if(exc==0){
					p = 1/(1+exp(-0.5*(chi+l0+log(sigma))));
					e_b = p*mu;
					e_b2 = p*(pow(mu,2)+sigma);
					//Rprintf("p: %g, e_b: %g, e_b2: %g\n",p,e_b,e_b2);
				}else{
					p = 0;
					e_b = mu;
					e_b2 = pow(mu,2);
					//Rprintf("p: %g, e_b: %g, e_b2: %g\n",p,e_b,e_b2);					
				}

				me(model,i,j)->p_sums = me(model,i,j)->p_sums + p;
				if(p>1-1e-10){
					me(model,i,j)->entropy = me(model,i,j)->entropy - p*log(p) + (1-p) + 0.5*p*log(2*3.14159*sigma);
				}else if(p<1e-10){
					me(model,i,j)->entropy = me(model,i,j)->entropy + p - (1-p)*log(1-p) + 0.5*p*log(2*3.14159*sigma);					
				} else {
					me(model,i,j)->entropy = me(model,i,j)->entropy - p*log(p) - (1-p)*log(1-p) + 0.5*p*log(2*3.14159*sigma);
				}
				me(model,i,j)->v_sums_correct = me(model,i,j)->v_sums_correct + (pow(e_b,2)-e_b2)*(model->data.x_sum_sq[k]);
					
				daxpy_w(model->data.n,xc(model,k),me(model,i,j)->resid_vec,me(model,i,j)->e_beta[k]-e_b);

				me(model,i,j)->beta_mu[k] = mu;
				me(model,i,j)->beta_sigma[k] = sigma;
				me(model,i,j)->beta_chi[k] = mu/sqrt(sigma);
				me(model,i,j)->e_beta[k] = e_b;
				me(model,i,j)->e_beta_sq[k] = e_b2;
				me(model,i,j)->beta_p[k] = p;
			}
			break;
		case LOGISTIC:
			//run logistic updates
			for(l=0;l< model->data.m ;l++){
				k = (&(model->data.ordering[i]))->col[l];
				for(t=0; t< model->data.n;t++){
					me(model,i,j)->x_w[t] = ((xc(model,k))[t])*(me(model,i,j)->w_vec[t]);
				}

				//k = l;
				//Rprintf("k: %d\n",k);
				exc = model->control_param.exclude[k];
				//exc = 1;
				//Rprintf("exc: %d\n",exc);
				//sigma = 1/((1/me(model,i,j)->sigma_e)*(model->data.x_sum_sq[k]));
				ddot_w(model->data.n,me(model,i,j)->x_w,xc(model,k),&prec);
				sigma = 1/prec;
				//Rprintf("sigma: %g\n",sigma);
				ddot_w(model->data.n,me(model,i,j)->x_w,me(model,i,j)->resid_vec,&mu);
				//Rprintf("mu: %g\n",mu);
				mu = mu + prec*(me(model,i,j)->e_beta[k]);
				//Rprintf("mu: %g\n",mu);
				mu = mu/prec;
				//Rprintf("mu: %g\n",mu);

				chi = pow(mu,2)/sigma;
				//Rprintf("chi: %g, e_beta[%d]: %g\n",chi,k,me(model,i,j)->e_beta[k]);
				if(exc==0){
					p = 1/(1+exp(-0.5*(chi+l0+log(sigma))));
					e_b = p*mu;
					e_b2 = p*(pow(mu,2)+sigma);
					//Rprintf("p: %g, e_b: %g, e_b2: %g\n",p,e_b,e_b2);
				}else{
					p = 0;
					e_b = mu;
					e_b2 = pow(mu,2);
					//Rprintf("p: %g, e_b: %g, e_b2: %g\n",p,e_b,e_b2);					
				}

				me(model,i,j)->p_sums = me(model,i,j)->p_sums + p;
				if(p>1-1e-10){
					me(model,i,j)->entropy = me(model,i,j)->entropy - p*log(p) + (1-p) + 0.5*p*log(2*3.14159*sigma);
				}else if(p<1e-10){
					me(model,i,j)->entropy = me(model,i,j)->entropy + p - (1-p)*log(1-p) + 0.5*p*log(2*3.14159*sigma);					
				} else {
					me(model,i,j)->entropy = me(model,i,j)->entropy - p*log(p) - (1-p)*log(1-p) + 0.5*p*log(2*3.14159*sigma);
				}
				//me(model,i,j)->v_sums_correct = me(model,i,j)->v_sums_correct + (pow(e_b,2)-e_b2)*(model->data.x_sum_sq[k]);
					
				daxpy_w(model->data.n,xc(model,k),me(model,i,j)->resid_vec,me(model,i,j)->e_beta[k]-e_b);
				daxpy_w(model->data.n,xc(model,k),me(model,i,j)->pred_vec_new,e_b - me(model,i,j)->e_beta[k]);

				me(model,i,j)->beta_mu[k] = mu;
				me(model,i,j)->beta_sigma[k] = sigma;
				me(model,i,j)->beta_chi[k] = mu/sqrt(sigma);
				me(model,i,j)->e_beta[k] = e_b;
				me(model,i,j)->e_beta_sq[k] = e_b2;
				me(model,i,j)->beta_p[k] = p;
			}

			break;

	}


}



void update_beta_marg(struct model_marg_struct * model, int * use_vec,int cv){

	int k,l,exc,t;
	double mu, sigma,prec, p, e_b;
  //double l0, e_b2, chi;
  
	//l0 = model->control_param.l0_path[j];
	//Rprintf("l0: %g\n",l0);
	//Rprintf("m: %d\n",model->data.m);

	switch (model->control_param.regressType){


		case LINEAR:
			//run linear updates
			for(l=0;l<cv;l++){
				//k = (&(model->data.ordering[i]))->col[l];
				k = use_vec[l];
				//k = l;
				//Rprintf("k: %d\n",k);
				exc = model->control_param.exclude[k];
				//Rprintf("exc: %d\n",exc);
				ddot_w(model->data.n,xcm(model,k),model->model_param.resid_vec,&mu);
				//Rprintf("mu: %g\n",mu);
				mu = mu + (model->data.x_sum_sq[k])*(model->model_param.beta_mu[k]);
				//Rprintf("mu: %g\n",mu);
				mu = mu/model->data.x_sum_sq[k];
				//Rprintf("mu: %g\n",mu);
				sigma = 1/((1/model->model_param.sigma_e)*(model->data.x_sum_sq[k]));
				//Rprintf("sigma: %g\n",sigma);
				//chi = pow(mu,2)/sigma;
				//Rprintf("chi: %g, e_beta[%d]: %g\n",chi,k,me(model,i,j)->e_beta[k]);
				if(exc==0){
					//p = 1/(1+exp(-0.5*(chi+l0+log(sigma))));
					p =0;
					//e_b = p*mu;
					e_b = mu;
					//e_b2 = p*(pow(mu,2)+sigma);
					//e_b2 = pow(mu,2);
					//Rprintf("p: %g, e_b: %g, e_b2: %g, k: %d\n",p,e_b,e_b2,k);
				}else{
					p = 0;
					e_b = mu;
					//e_b2 = pow(mu,2);
					//Rprintf("p: %g, e_b: %g, e_b2: %g, k: %d\n",p,e_b,e_b2,k);					
				}

				//me(model,i,j)->p_sums = me(model,i,j)->p_sums + p;
				//me(model,i,j)->v_sums_correct = me(model,i,j)->v_sums_correct + (pow(e_b,2)-e_b2)*(model->data.x_sum_sq[k]);
					
				daxpy_w(model->data.n,xcm(model,k),model->model_param.resid_vec,model->model_param.beta_mu[k]-e_b);
        //Rprintf("daxpy\n");
				model->model_param.beta_mu[k] = mu;
        //Rprintf("mu\n");
				model->model_param.beta_sigma[k] = sigma;
        //Rprintf("sigma\n");
				model->model_param.beta_chi[k] = mu/sqrt(sigma);
        //Rprintf("chi\n");
				model->model_param.beta_p[k] = p;
        //Rprintf("p\n");
			}
			break;
		case LOGISTIC:
			//run logistic updates
			for(l=0;l<cv;l++){
				//k = (&(model->data.ordering[i]))->col[l];
				k = use_vec[l];
				for(t=0; t< model->data.n;t++){
					model->model_param.x_w[t] = ((xcm(model,k))[t])*(model->model_param.w_vec[t]);
				}

				//k = l;
				//Rprintf("k: %d\n",k);
				exc = model->control_param.exclude[k];
				//exc = 1;
				//Rprintf("exc: %d\n",exc);
				//sigma = 1/((1/me(model,i,j)->sigma_e)*(model->data.x_sum_sq[k]));
				ddot_w(model->data.n,model->model_param.x_w,xcm(model,k),&prec);
				sigma = 1/prec;
				//Rprintf("sigma: %g\n",sigma);
				ddot_w(model->data.n,model->model_param.x_w,model->model_param.resid_vec,&mu);
				//Rprintf("mu: %g\n",mu);
				mu = mu + prec*(model->model_param.beta_mu[k]);
				//Rprintf("mu: %g\n",mu);
				mu = mu/prec;
				//Rprintf("mu: %g\n",mu);

				//chi = pow(mu,2)/sigma;
				//Rprintf("chi: %g, e_beta[%d]: %g\n",chi,k,me(model,i,j)->e_beta[k]);
				if(exc==0){
					//p = 1/(1+exp(-0.5*(chi+l0+log(sigma))));
					p =0;
					//e_b = p*mu;
					e_b = mu;
					//e_b2 = p*(pow(mu,2)+sigma);
					//e_b2 = pow(mu,2);
					//Rprintf("p: %g, e_b: %g, e_b2: %g, sigma: %g, k: %d\n",p,e_b,e_b2,sigma,k);
				}else{
					p = 0;
					e_b = mu;
					//e_b2 = pow(mu,2);
					//Rprintf("p: %g, e_b: %g, e_b2: %g, sigma: %g, k: %d\n",p,e_b,e_b2,sigma,k);					
				}

				//me(model,i,j)->p_sums = me(model,i,j)->p_sums + p;
				//me(model,i,j)->v_sums_correct = me(model,i,j)->v_sums_correct + (pow(e_b,2)-e_b2)*(model->data.x_sum_sq[k]);
					
				daxpy_w(model->data.n,xcm(model,k),model->model_param.resid_vec,model->model_param.beta_mu[k]-e_b);
				daxpy_w(model->data.n,xcm(model,k),model->model_param.pred_vec_new,e_b - model->model_param.beta_mu[k]);

				model->model_param.beta_mu[k] = mu;
				model->model_param.beta_sigma[k] = sigma;
				model->model_param.beta_chi[k] = mu/sqrt(sigma);
				model->model_param.beta_p[k] = p;
			}

			break;

	}


}


void update_error(struct model_struct * model, int i, int j){

	int t;
	double U;
	double nd = (double) model->data.n;
	
	switch(model->control_param.regressType){

		case LINEAR:

			ddot_w(model->data.n,me(model,i,j)->resid_vec,me(model,i,j)->resid_vec,&U);
			U = U - me(model,i,j)->v_sums_correct;
			U = U/nd;
			me(model,i,j)->sigma_e = U;
      //me(model,i,j)->sigma_e = 1.0;
                        //Rprintf("no segfault\n");
			if(!R_FINITE(U)){
				free_model(model);
				//Rprintf("segfault\n");
				error("Penalized linear solution does not exist.\n");
				//error("uh oh\n");
			}

			break;
		
		case LOGISTIC:
			////
			for(t=0;t<model->data.n;t++){
				
				me(model,i,j)->mu_vec[t] = 1/(1+exp(-me(model,i,j)->pred_vec_new[t]));
				me(model,i,j)->w_vec[t] = me(model,i,j)->mu_vec[t]*(1-me(model,i,j)->mu_vec[t]);
				me(model,i,j)->resid_vec[t] = (model->data.y[t]-me(model,i,j)->mu_vec[t])/me(model,i,j)->w_vec[t];
				me(model,i,j)->pred_vec_old[t] = me(model,i,j)->pred_vec_new[t];
				if(me(model,i,j)->mu_vec[t]==1 || me(model,i,j)->mu_vec[t] ==0){
					//Rprintf("OVERFIT\n");
					free_model(model);
					error("Penalized logistic solution does not exist.\n");
				}

			}

			break;

	}



}


void update_error_marg(struct model_marg_struct * model){

	int t;
	double U;
	double nd = (double) model->data.n;
	
	switch(model->control_param.regressType){

		case LINEAR:

			ddot_w(model->data.n,model->model_param.resid_vec,model->model_param.resid_vec,&U);
			//U = U - me(model,i,j)->v_sums_correct;
			U = U/nd;
			model->model_param.sigma_e = U;
			if(!R_FINITE(U)){error("Penalized linear solution does not exist.\n");}

			break;
		
		case LOGISTIC:
			////
			for(t=0;t<model->data.n;t++){
				
				model->model_param.mu_vec[t] = 1/(1+exp(-model->model_param.pred_vec_new[t]));
				model->model_param.w_vec[t] = model->model_param.mu_vec[t]*(1-model->model_param.mu_vec[t]);
				model->model_param.resid_vec[t] = (model->data.y[t]-model->model_param.mu_vec[t])/model->model_param.w_vec[t];
				model->model_param.pred_vec_old[t] = model->model_param.pred_vec_new[t];
				if(model->model_param.mu_vec[t]==1 || model->model_param.mu_vec[t] ==0){
					//Rprintf("OVERFIT\n");
					error("Penalized logistic solution does not exist.\n");
				}

			}

			break;

	}



}


//void update_p_beta(struct model_struct * model, int i, int j){

//	double md = (double) model->data.m;
//	double pd = (double) model->data.p;
//	double p_beta, l0;
//
//	p_beta = (me(model,i,j)->p_sums)/(md-pd);
//	//Rprintf("pd:%g, p_sums: %g, p_beta: %g, md: %g\n",pd,me(model,i,j)->p_sums,p_beta,md);
//	l0 = 2*(log(p_beta)-log(1-p_beta));
//	me(model,i,j)->p_max = p_beta;
//	me(model,i,j)->l0_max = l0;
//	
//	
//}

void update_lb(struct model_struct * model, int i, int j){



	double lba;
	double nd = (double) model->data.n;
	double md = (double) model->data.m;
	double p_beta;
	//if(model->control_param.max_pb==1){
	//	p_beta = me(model,i,j)->p_max;
	//}else{
		p_beta = model->control_param.pb_path[j];
	//}
	int t;

	switch(model->control_param.regressType){
		
		case LINEAR:

			lba = -0.5*nd*(log(2*3.14159*me(model,i,j)->sigma_e)+1);
			lba = lba + log(p_beta)*(me(model,i,j)->p_sums);
			lba = lba + log(1-p_beta)*(md - me(model,i,j)->p_sums);
			lba = lba + me(model,i,j)->entropy;
			//Rprintf("Entropy: %g\n",me(model,i,j)->entropy);
			me(model,i,j)->lb = lba;

			break;
		
		case LOGISTIC:
			////

			//lba = -0.5*(log(me(model,i,j)->sigma_e)+1);
			ddot_w(model->data.n,model->data.y,me(model,i,j)->pred_vec_new,&lba);
			for(t=0;t<model->data.n;t++){
				lba = lba + log(1-me(model,i,j)->mu_vec[t]);
			}
			lba = lba + log(p_beta)*(me(model,i,j)->p_sums);
			lba = lba + log(1-p_beta)*(md - me(model,i,j)->p_sums);
			lba = lba + me(model,i,j)->entropy;
			me(model,i,j)->lb = lba;

			break;

	}



}


void update_lb_marg(struct model_marg_struct * model){



	double lba;
	//double nd = (double) model->data.n;
	//double md = (double) model->data.m;
	int t;

	switch(model->control_param.regressType){
		
		case LINEAR:

			lba = -0.5*(log(model->model_param.sigma_e)+1);
			//lba = lba + (model->control_param.pb_path[j])*(me(model,i,j)->p_sums);
			//lba = lba + (1-model->control_param.pb_path[j])*(md - me(model,i,j)->p_sums);
			//me(model,i,j)->lb = lba;
			model->model_param.lb = lba;

			break;
		
		case LOGISTIC:
			////

			//lba = -0.5*(log(me(model,i,j)->sigma_e)+1);
			ddot_w(model->data.n,model->data.y,model->model_param.pred_vec_new,&lba);
			for(t=0;t<model->data.n;t++){
				lba = lba + log(1-model->model_param.mu_vec[t]);
			}
			//lba = lba + (model->control_param.pb_path[j])*(me(model,i,j)->p_sums);
			//lba = lba + (1-model->control_param.pb_path[j])*(md - me(model,i,j)->p_sums);
			//me(model,i,j)->lb = lba;
			model->model_param.lb = lba;

			break;

	}



}


void run_vbsr(struct model_struct * model){	
	int i,j;
	double tol=1;
	double lb_old;
	int count = 0;
	//#pragma omp parallel for private(i,j,count,tol,lb_old)
	for (i=0;i < model->control_param.n_orderings;i++){
		for(j=0;j < model->control_param.path_length;j++){
			if(j>0){
				//copy the previous path to the new path
				copy_model_state(model,i,j);
				//Rprintf("Copied model state...\n");
			}
			while(fabs(tol) > model->control_param.eps && count < model->control_param.maxit){
				
				me(model,i,j)->p_sums = 0;
				me(model,i,j)->v_sums_correct = 0;
				me(model,i,j)->entropy = 0;
				lb_old = me(model,i,j)->lb;
				//Rprintf("Updating beta...\n");
				update_beta(model,i,j);

				//if(model->control_param.max_pb==1){
				//	update_p_beta(model,i,j);
				//}
				//Rprintf("Updating error...\n");
				update_error(model,i,j);
				//Rprintf("Updating lower bound...\n");
				update_lb(model,i,j);
				tol = lb_old - me(model,i,j)->lb;
				count = count+1;

			}
			//Rprintf("lb: %g,i: %d, j: %d\n",lb_old,i,j);			
			if(count>=model->control_param.maxit){
				Rprintf("Maximum iterations exceeded!\n");
			}
			count =0;
			tol = 1;
		}
	}

}

void run_marg(struct model_marg_struct * model){	
	int j,k,cv,q;
	double tol=1;
	double lb_old;
	int count = 0;
	//#pragma omp parallel for
	//for (i=0;i < model->control_param.n_orderings;i++){
	//	for(j=0;j < model->control_param.path_length;j++){
	//		if(j>0){
				//copy the previous path to the new path
	//			copy_model_state(model,i,j);
				//Rprintf("Copied model state...\n");
	//		}

	for (j=0;j<model->data.m;j++){
		if(model->control_param.exclude[j]==1){count = count +1;}
	}
	cv = count+1;
	count = 0;
	int * use_vec = (int *) malloc(sizeof(int)*(cv));
	for(j=0;j<model->data.m;j++){
		if(model->control_param.exclude[j]==1){
			use_vec[count+1] = j;
			count = count +1;
		}
	}
	

	for(j=0;j<model->data.m;j++){
		
		if(model->control_param.exclude[j]==0){
			for(k=0;k<model->data.n;k++){
				switch(model->control_param.regressType){
					case LINEAR:
						model->model_param.resid_vec[k] = model->data.y[k];
						
						break;
					case LOGISTIC:
						model->model_param.resid_vec[k] = (model->data.y[k] - 0.5)/(0.25);
						model->model_param.w_vec[k] = 0.25;
						model->model_param.mu_vec[k] = 0.5;
						model->model_param.pred_vec_old[k] = 0;
						model->model_param.pred_vec_new[k] = 0;
						break;
				}

			}
			tol =1;
			count = 0;
			//Rprintf("JJJJJ: %d\n",j);
			use_vec[0] = j;
			//set fixed covariates to zero!
			for (q=1;q<cv;q++){
				model->model_param.beta_mu[use_vec[q]] = 0;
				model->model_param.beta_sigma[use_vec[q]] = 0;
				model->model_param.beta_chi[use_vec[q]] = 0;
				model->model_param.beta_p[use_vec[q]] = 0;
			}
			model->model_param.lb = 0;	
			while(fabs(tol) > model->control_param.eps && count < model->control_param.maxit){
				//me(model,i,j)->p_sums = 0;
				//me(model,i,j)->v_sums_correct = 0;
				lb_old = model->model_param.lb;
				//Rprintf("Updating beta...\n");
				update_beta_marg(model,use_vec,cv);
				//Rprintf("Updating error...\n");
				update_error_marg(model);
				//Rprintf("Updating lower bound...\n");
				update_lb_marg(model);
				tol = lb_old - model->model_param.lb;
				count = count+1;
			}
			if(count>=model->control_param.maxit){
				warning("Maximum iterations exceeded!\n");
			}
		}
	}
	//Rprintf("lb: %g,i: %d, j: %d\n",lb_old,i,j);			

	free(use_vec);
	//count =0;
	//tol = 1;
	
}


void identify_unique(double * lb_t, double * post_p, int n,double tol){
	int i,j,count;
	double tv;
	count =0;

	for(i=0;i<n-1;i++){
		for(j=i+1;j<n;j++){
			if(i!=j){
				tv = fabs(lb_t[i]-lb_t[j]);
				if(tv < tol){
					post_p[j] = 0;
				}
			}
		}
	}

	
	tv =0;
	for(i=0;i<n;i++){
		if(post_p[i]>0){count=count+1;}
		tv = post_p[i]+tv;
	}
	for(i=0;i<n;i++){
		post_p[i]=post_p[i]/tv;
		//if(post_p[i]>0){Rprintf("post_prob[%d]: %g\t",i,post_p[i]);}
	}
	//Rprintf("Identified: %d unique models\n",count);
}

void compute_bma_correct(struct model_struct * model,int k,double * post_prob,double * s_bma,int j){
	int t,l;
	double corv;
	s_bma[0] = 0;

	//t ord ind
	//l ord ind
	//k marker ind
	//j path ind
	for (t=0;t<model->control_param.n_orderings;t++){
		if(post_prob[t] > 0){
			s_bma[0] = s_bma[0] + pow(post_prob[t],2);
		}
	}


	for(t=0;t<model->control_param.n_orderings-1;t++){
		for(l=t+1;l<model->control_param.n_orderings;l++){
			if(post_prob[t]>0 && post_prob[l]>0){
				daxpy_w(model->data.n,xc(model,k),me(model,t,j)->resid_vec,me(model,t,j)->e_beta[k]);
				daxpy_w(model->data.n,xc(model,k),me(model,l,j)->resid_vec,me(model,l,j)->e_beta[k]);
				cor(me(model,t,j)->resid_vec, me(model,l,j)->resid_vec, model->data.one_vec,&corv,model->data.n);
				daxpy_w(model->data.n,xc(model,k),me(model,t,j)->resid_vec,-me(model,t,j)->e_beta[k]);
				daxpy_w(model->data.n,xc(model,k),me(model,l,j)->resid_vec,-me(model,l,j)->e_beta[k]);
				s_bma[0] = s_bma[0] + 2*post_prob[t]*post_prob[l]*(corv);
				//if(j==2 && k==0){Rprintf("correction: %g %g %g %g\n",s_bma[0],corv,post_prob[t],post_prob[l]);}
			}
		}
	}
	//Rprintf("correction: %g\n",s_bma[0]);




}


void collapse_results(struct model_struct * model,
						double * beta_chi_mat, 
						double * beta_mu_mat, 
						double * beta_sigma_mat, 
						double * e_beta_mat, 
						double * beta_p_mat, 
						double * lb_mat, 
						double * kl_mat){
	
	int i,j,k;
	double max_v,bc,bm,bs,eb,bp,Z,s_bma;
	double * post_prob = (double *) malloc(sizeof(double)*model->control_param.n_orderings);
	double * lb_t = (double *) malloc(sizeof(double)*model->control_param.n_orderings);
	//max_v = -1e100;
	int w_max;
	//if(model->control_param.max_pb==1){
	//	for(i=0;i<model->control_param.n_orderings;i++){
	//		p_est[i]=me(model,i,0)->p_max;
	//	}
	//}
	
	switch(model->control_param.estType){
		
		
		case BMA:
			for(j=0;j<model->control_param.path_length;j++){	
				max_v = me(model,0,j)->lb;
				w_max = 0;
				Z =0;
				for(i=0;i<model->control_param.n_orderings;i++){
					if(me(model,i,j)->lb > max_v){
						max_v = me(model,i,j)->lb;
						w_max = i;
					}
					lb_mat[(model->control_param.n_orderings)*(j)+i] = me(model,i,j)->lb;
					lb_t[i] = me(model,i,j)->lb;
				}
				for(i=0;i<model->control_param.n_orderings;i++){
					Z = Z + exp(me(model,i,j)->lb-max_v);
				}
				for(i=0;i<model->control_param.n_orderings;i++){
					post_prob[i] = exp(me(model,i,j)->lb-max_v)/Z;
					//Rprintf("post_prob[%d]: %g\t",i,post_prob[i]);
					//lb_mat[(model->control_param.n_orderings)*(j)+i] = post_prob[i];
					//bm = bm + post_prob*me(model,i,j)->beta_chi
				}
				//Rprintf("\n");

				identify_unique(lb_t,post_prob,model->control_param.n_orderings,model->control_param.eps*10);
				
				for(k=0;k<model->data.m;k++){
					bc =0;
					bm =0;
					bs=0;
					eb=0;
					bp=0;
					switch(model->control_param.errType){
					case APPR:
						s_bma = 1;
						break;
					case EXACT:
						compute_bma_correct(model,k,post_prob,&s_bma,j);
						break;
					default:
						Rprintf("BMA computation not specified!\n");
						break;
					}

					for(i=0;i<model->control_param.n_orderings;i++){
						bc = bc+ post_prob[i]*me(model,i,j)->beta_chi[k];
						bm = bm+ post_prob[i]*me(model,i,j)->beta_mu[k];
						bs = bs+ post_prob[i]*me(model,i,j)->beta_sigma[k];
						eb = eb+ post_prob[i]*me(model,i,j)->e_beta[k];
						bp = bp+ post_prob[i]*me(model,i,j)->beta_p[k];
					}
					beta_chi_mat[(model->data.m)*(j)+k] = bc/sqrt(s_bma);
					beta_mu_mat[(model->data.m)*(j)+k] = bm;
					beta_sigma_mat[(model->data.m)*(j)+k] = bs;
					e_beta_mat[(model->data.m)*(j)+k] = eb;
					beta_p_mat[(model->data.m)*(j)+k] = bp;
				}				
			}
					
			break;
			
		case MAXIMAL:
			////
			for(j=0;j<model->control_param.path_length;j++){	
				max_v = me(model,0,j)->lb;
				w_max = 0;
				for(i=0;i<model->control_param.n_orderings;i++){
					if(me(model,i,j)->lb > max_v){
						max_v = me(model,i,j)->lb;
						w_max = i;
					}
					lb_mat[(model->control_param.n_orderings)*(j)+i] = me(model,i,j)->lb;
				}
				for(k=0;k<model->data.m;k++){
					beta_chi_mat[(model->data.m)*(j)+k] = me(model,w_max,j)->beta_chi[k];
					beta_mu_mat[(model->data.m)*(j)+k] = me(model,w_max,j)->beta_mu[k];
					beta_sigma_mat[(model->data.m)*(j)+k] = me(model,w_max,j)->beta_sigma[k];
					e_beta_mat[(model->data.m)*(j)+k] = me(model,w_max,j)->e_beta[k];
					beta_p_mat[(model->data.m)*(j)+k] = me(model,w_max,j)->beta_p[k];			
				}				
			}
			break;
			
		default:
			////
			break;
			
			
	}
	free(post_prob);
	free(lb_t);

}

void collapse_results_marg(struct model_marg_struct * model,
							double * beta_chi,
							double * beta_mu,
							double * beta_sigma,
							double * beta_p,
							double * lb){
							
							
	int j;
	*lb = model->model_param.lb;
	for(j=0;j<model->data.m;j++){
		beta_chi[j]=model->model_param.beta_chi[j];
		beta_mu[j]=model->model_param.beta_mu[j];
		beta_sigma[j]=model->model_param.beta_sigma[j];
		beta_p[j]=model->model_param.beta_p[j];
	}
}

void run_marg_analysis(double * eps,
			int * exclude,
			int * maxit,
			int * regress,
			int * scale,
			double * X,
			double * y,
			double * var_y,
			int * n,
			int * m,
			double * beta_chi,
			double * beta_mu,
			double * beta_sigma,
			double * beta_p,
			double * lb){
			
			
	struct model_marg_struct model;
	//Rprintf("Initializing marginal analysis...\n");
	initialize_model_marg(eps,exclude,maxit,regress,scale,X,y,var_y,n,m,&model);
	//Rprintf("Initialized marginal model...\n");
	run_marg(&model);
	//Rprintf("Model run...\n");
	collapse_results_marg(&model,beta_chi,beta_mu,beta_sigma,beta_p,lb);
	free_model_marg(&model);
			
			
}

void run_vbsr_wrapper(double * eps, 
			double * l0_path, 
			double * pb_path, 
			int * exclude, 
			double * penalty_factor, 
			int * maxit, 
			int * path_length, 
			int * n_orderings, 
			int * regress, 
			int * scale, 
			int * est, 
			int * error, 
			double * kl,
			int * approx,
			int * total_replicates,
			double * X, 
			double * y, 
			double * var_y, 
			int * n, 
			int * m, 
			int * ordering_mat, 
			double * beta_chi_mat, 
			double * beta_mu_mat, 
			double * beta_sigma_mat, 
			double * e_beta_mat, 
			double * beta_p_mat, 
			double * lb_mat, 
			double * kl_mat,
			int * nthreads){
						
	
	struct model_struct model;
	//omp_set_num_threads(*nthreads);
	//Rprintf("nthreads: %d, nthreads_o: %d\n",*nthreads,omp_get_max_threads());
	//Rprintf("Initializing model...\n");
	initialize_model(eps,l0_path,pb_path,exclude,penalty_factor,maxit,path_length,n_orderings,regress,scale,est,error,kl,approx,total_replicates,X, y, var_y, n, m,ordering_mat,&model);
	//Rprintf("Initialized model...\n");
	run_vbsr(&model);
	//Rprintf("Model run...\n");
	collapse_results(&model,beta_chi_mat, beta_mu_mat, beta_sigma_mat, e_beta_mat, beta_p_mat, lb_mat, kl_mat);
	//Rprintf("Results computed..\n");
	free_model(&model);

}




