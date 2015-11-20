#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string.h>


double likelihood(double *y_obs, double* y_model);
void my_model(double *y_init, double *x_obs, double a, double b, double c, double d);
int minimo (double *arr, int nsteps, int nburn);
double std_dev(double *arr, int n_steps, int n_burn);

int main(int argc, char **argv){
 
  
  

  int n_steps = atoi(argv[1]);
  int n_burn = atoi(argv[2]);
  int npoints = n_steps;

  double *x_obs;
  double *y_obs;
  double *a_walk;
  double *b_walk;
  double *c_walk;
  double *d_walk;
  double *l_walk;
  double *y_init;
  double *y_prime;
  FILE *in;
  double var0;
  double var1;
  double var2;
  double var3;
  double var4;

  x_obs  = malloc(npoints * sizeof(double));
  y_obs  = malloc(npoints * sizeof(double));
  a_walk = malloc(npoints * sizeof(double));
  b_walk = malloc(npoints * sizeof(double));
  c_walk = malloc(npoints * sizeof(double));
  d_walk = malloc(npoints * sizeof(double));
  l_walk = malloc(npoints * sizeof(double));
  y_init = malloc(npoints * sizeof(double));
  y_prime = malloc(npoints * sizeof(double));

  in = fopen("monthrg.dat", "r");
  if(!in){
    printf("problems opening the file %s\n", "monthrg.dat");
    exit(1);
  }
  
  int i;
  for(i=0;i<4141;i++){
    fscanf(in,"%lf %lf %lf %lf %lf\n", &var0, &var1, &var2, &var3, &var4);
    if (var3!=-99){
      x_obs[i] = var0 + var1/12.0;
      y_obs[i] = var3;
      // printf("%f %f\n",x_obs[i], var3);
    }
      
   
  }
  fclose(in);
  
  char output1[512];
  char output2[512];
  char* line = "_";
  char* mcmc_solar = "mcmc_solar";
  char* dat = ".dat";
  char name[512] = "";
  strcat(name,mcmc_solar);
  strcat(name,line);
  strcat(name,argv[1]);
  strcat(name,line);
  strcat(name,argv[2]);
  strcat(name,dat);
  //printf("%s\n", name);

  srand((unsigned)time(NULL));
		   

  
  FILE *fp;
  fp = fopen(name, "ab+");
  
 
  
  
  
  
  
  a_walk[0]=((double)rand() / (double) (RAND_MAX/50.0));
  b_walk[0]=((double)rand() / (double) (RAND_MAX/50.0));
  c_walk[0]=((double)rand() / (double) (RAND_MAX/50.0));
  d_walk[0]=((double)rand() / (double) (RAND_MAX/50.0));
  
  // printf("%f %f %f %f\n",a_walk[0], b_walk[0], c_walk[0], d_walk[0]);
  my_model (y_init, x_obs, a_walk[0], b_walk[0], c_walk[0], d_walk[0]);
  l_walk[0] = likelihood(y_obs, y_init);
  // printf("%f\n",l_walk[0]);
  

  int j;
  
  
   const gsl_rng_type * T;
   gsl_rng * r;
   gsl_rng_env_setup();
   T = gsl_rng_default;
   r = gsl_rng_alloc(T);

   double a_prime;
   double b_prime;
   double c_prime;
   double d_prime;
   double l_prime;
   double l_init;
   double gamma;
   double beta;
  
   for (j = 0; j<n_steps; j++){

    a_prime =(double)gsl_ran_gaussian(r, 0.1)+a_walk[j];
    b_prime =(double)gsl_ran_gaussian(r, 0.1)+b_walk[j];
    c_prime =(double)gsl_ran_gaussian(r, 0.1)+c_walk[j];
    d_prime =(double)gsl_ran_gaussian(r, 0.1)+d_walk[j];
    
    // printf("%f %f %f %f\n",a_walk[j], b_walk[j], c_walk[j], d_walk[j]);

    my_model(y_init, x_obs, a_walk[j], b_walk[j], c_walk[j], d_walk[j]);
    my_model(y_prime, x_obs, a_prime, b_prime, c_prime, d_prime);

    l_prime = likelihood(y_obs, y_prime);
    l_init = likelihood(y_obs, y_init);

    gamma = -l_prime + l_init;
    // printf("%f\n",gamma);

    if(gamma>=0.0){

      a_walk[j+1] = a_prime;
      b_walk[j+1] = b_prime;
      c_walk[j+1] = c_prime;
      d_walk[j+1] = d_prime;
      l_walk[j+1] = l_prime;

    }
    else{
      
      beta = ((double)rand() / (double) (RAND_MAX/50.0));
      
      if(beta<=exp(gamma)){

      a_walk[j+1] = a_prime;
      b_walk[j+1] = b_prime;
      c_walk[j+1] = c_prime;
      d_walk[j+1] = d_prime;
      l_walk[j+1] = l_prime;

      }
      else{

      a_walk[j+1] = a_walk[j];
      b_walk[j+1] = b_walk[j];
      c_walk[j+1] = c_walk[j];
      d_walk[j+1] = d_walk[j];
      l_walk[j+1] = l_init;

      }
      

      
      
    }
    // printf("%f %f %f %f\n",a_walk[j], b_walk[j], c_walk[j], d_walk[j]);
    
  }
  double *a_burned;
  double *b_burned;
  double *c_burned;
  double *d_burned;
  double *l_burned;
  double *x_burned;

  a_burned = malloc(npoints * sizeof(double));
  b_burned = malloc(npoints * sizeof(double));
  c_burned = malloc(npoints * sizeof(double));
  d_burned = malloc(npoints * sizeof(double));
  l_burned = malloc(npoints * sizeof(double));
  
  
  int co;

  for (co = 0; co < n_steps-n_burn; co++){

    a_burned[co] = a_walk[n_burn + co];
    b_burned[co] = b_walk[n_burn + co];
    c_burned[co] = c_walk[n_burn + co];
    d_burned[co] = d_walk[n_burn + co];
    l_burned[co] = l_walk[n_burn + co];
    fprintf(fp, "%f %f %f %f %f\n",a_burned[co], b_burned[co], c_burned[co], d_burned[co], l_burned[co]);
    //printf("%f %f %f %f %f\n",a_walk[co], b_walk[co], c_walk[co], d_walk[co], l_walk[co]);
   
  }

  int max_likelihood_id = minimo(l_burned, n_steps, n_burn);
  // printf("%i\n", max_likelihood_id);
  double best_a = a_burned[max_likelihood_id];
  double best_b = b_burned[max_likelihood_id];
  double best_c = c_burned[max_likelihood_id];
  double best_d = d_burned[max_likelihood_id];
  // printf("%f %f %f %f\n",a_burned[max_likelihood_id], b_burned[max_likelihood_id], c_burned[max_likelihood_id], d_burned[max_likelihood_id]);

  printf("a = %f +/- %f\n",best_a, std_dev(a_burned, n_steps, n_burn));
  printf("b = %f +/- %f\n",best_b, std_dev(b_burned, n_steps, n_burn));
  printf("c = %f +/- %f\n",best_c, std_dev(c_burned, n_steps, n_burn));
  printf("d = %f +/- %f\n",best_d, std_dev(d_burned, n_steps, n_burn));

 
   
 fclose(fp);




  gsl_rng_free(r);
  return 0;
}





















double likelihood(double *y_obs, double* y_model){
  
  double chi_squared = 0.0;
  int counter;
 
  for (counter = 0; counter < 4141; counter++){
    chi_squared = chi_squared + pow((y_obs[counter]-y_model[counter]),2);
  }
  return (1.0/2.0)*chi_squared;
}

void my_model(double *y_init, double *x_obs, double a, double b, double c, double d){
  
  int counter2;
  double pi = 3.14159265359;
  
  for (counter2 = 0; counter2<4141; counter2++){

    y_init[counter2] = a*cos((2*pi/d)*x_obs[counter2] + b) + c;

  }
}

int minimo (double *arr, int nsteps, int nburn){

    int c = 0;
    double min = arr[0];
    int minimo = 0;

    for (c = 1; c<nsteps-nburn; c++){

      if (arr[c] < min){

	min = arr[c];
	minimo = c;

      
      }

    }
    return minimo;

}
  
double std_dev(double *arr, int n_steps, int n_burn){

  double mean = 0.0, sum_dev = 0.0;

  int k;
  int n = n_steps-n_burn;

  for (k=0; k<n; k++){
    
    mean+=arr[k];
    
  }
  mean = mean/n;
  for (k = 0; k<n; k++){

    sum_dev+=(arr[k]-mean)*(arr[k]-mean);
    
  }
  return pow(sum_dev/n,0.5);
  

}
  


  
  



  
  
 
    
 
    
        
     

