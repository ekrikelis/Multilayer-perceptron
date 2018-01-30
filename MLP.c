// implementation of a multilayer perceptron in C
// https://en.wikipedia.org/wiki/Multilayer_perceptron

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define d	2		//input dimension (arithmos eisodwn MLP)
#define K	3		//output dimension (arithmos katigoriwn (eksodwn) MLP)
#define H1	7		//number of nodes in the 1st hidden layer (arithmos krimenwn nevronwn gia 1o epipedo)
#define H2	2		//number of nodes in the 2nd hidden layer (arithmos krimenwn nevronwn gia 2o epipedo)

#define e	0.0001		//error rate ( minimum metavoli sfalmatos )
#define N	1000		//Number of training set ( arithmos protypou synolou ekpaideushs)
#define n	0.01		//descent step (vima kathodou)
#define MaxEpoxh	30000		//number of epoches (max)
#define T	1000		//number of test set (arithmos synolou elegxou)



//int flag=1;			//flag=1 logistic activation function (logistiki synartisi energopoihshs)
				//flag=2 linear activation function (grammiki synartisi energopoihshs)

double E;			//training error (synoliko sfalama ekpaideusis MLP)
double Epre;			//training error per epoche (sfalma ekpaideusis prohgoumenhs epoxis)



//initilization of MLP architecture
/*Pinakes polosewn*/
double b1[H1];		//array of bias of nodes of the 1st hidden layer - pinakas polosewn 1ou epipedou
double b2[H2];		//array of bias of nodes of the 2nd hidden layer - pinakas polosewn 2ou epipedou
double b3[K];		//array of bias of nodes of the ouput layer - pinakas polosewn eksodou

/*arrays of weights*/
double w1[H1][d];	//array of weights pinakas from input to  1st hidden layer
double w2[H2][H1];	//array of weights pinakas from 1st hidden to 2nd hidden layer
double w3[K][H2];	//array of weights pinakas from 2nd hidden to output

//Outputs - eksodoi
double y[K]; //output of MLP (eksodoi tou MLP)
double z1[H1];//output of 1st layer (eksodoi nevrwnwn 1ou epipedou)
double z2[H2];//output of the 2nd layer (eksodoi nevrwnwn 2ou epipedou)

//Error - sfalmata
double d1[H1];	//errors of the 1st layer (sfalmata nevrwnwn 1ou epipedou)
double d2[H2];	//errors of the 2nd layer (sfalmata nevrwnwn 2ou epipedou)
double d3[K];	//errors of the output	  (sfalmata nevrwnwn eksodou)

// the sum of producers of weights- atroisma paragwgwn varwn
double dw1[H1][d];	//eisodou->1ou epipedou 
double dw2[H2][H1];	//1ou-->2ou epipedou
double dw3[K][H2];	//2ou-->eksodou

//the sum of bias - athroisma paragwgwn polosewn
double db1[H1];		//nodes of 1st layer
double db2[H2];		//nodes of 2nd layer
double db3[K];		//nodes of output

//eswterika ginomena s=x*w
double	s1[H1];	
double	s2[H2];	
double	s3[K];

//arrays for the training set - pinakes protypwn ekpaideusis 
double x_train[N][d];		//input 
double t_train[N][K];		//target 

//arrays for the test set - pinakes protypwn elegxou 
double x_test[T][d];		//input 
double t_test[T][K];		//target

FILE* f1;
FILE* f2;
FILE* f3;

//----------synartiseis-----------------------------

//arxikopoihsh paragogwn polosewn sto 0
void initialiaze_parag_polosewn()
{
	int i;

	//paragwgoi polosewn 1ou epipedou
	for(i=0; i<H1; i++)
		db1[i]=0;

	//paragwgoi polosewn 2ou epipedou
	for(i=0; i<H2; i++)
		db2[i]=0;

	//paragwgoi polosewn eksodou
	for(i=0; i<K; i++)
		db3[i]=0;
	
		
}

//arxikopoihsh paragogwn varwn sto 0
void initialiaze_parag_varwn()
{
	//paragwgoi polosewn 1ou epipedou
	int i,j;
	for(i=0; i<H1; i++)
		for(j=0; j<d; j++)
			dw1[i][j]=0;


	//paragwgoi polosewn 2ou epipedou
	for(i=0; i<H2; i++)
		for(j=0; j<H1; j++)
			dw2[i][j]=0;

	//paragwgoi polosewn eksodou
	for(i=0; i<K; i++)
		for(j=0; j<H2; j++)
			dw3[i][j]=0;

	
		
}

int init=0;
// creation of a random number - dimiourgia tyxaiou arithmou
double GetRand()
{
 int i;
 while(init==0)
 {
    srand((unsigned)(time(0)));
    init=1;
 }
 
	  double rr = (   ((double)rand() / ((double)(RAND_MAX))*2-1.0));
	      
	  return(rr);
 
}






// initialization of bias - arxikopoihsh polosewn sto (-1,1)
void initialize_poloseis()
{
	int i;
	
	//poloseis 1ou epipedou
	for(i=0; i<H1; i++)
		b1[i]=GetRand();

	//poloseis 2ou epipedou
	for(i=0; i<H2; i++)
		b2[i]=GetRand();
	
	//poloseis eksodou
	for(i=0; i<K; i++)
		b3[i]=GetRand();
	
}

//initialization of weights in (-1,1)
void initialize_varwn()
{
	int i,j;

	//varoi eisodou-->1ou epipedou
	for(i=0; i<H1; i++)
		for(j=0; j<d; j++)
			w1[i][j]=GetRand();

	//varoi 1ou-->2ou epipedou
	for(i=0; i<H2; i++)
		for(j=0; j<H1; j++)
			w2[i][j]=GetRand();

	
	//varoi 2ou-->eksodou
	for(i=0; i<K; i++)
		for(j=0; j<H2; j++)
			w3[i][j]=GetRand();
	


}


// activation function
double f(double u)
{
	//if (flag==1)		//logistiki synartisi energopoihshs
		return (1/(1+exp(-u)));	

	//if(flag==2)		//grammiki synartisi energopoihshs
		return (u);
}


//forward_pass---euthy perasma


//prosoxi isws den xreiazetai o pinakas s1,s2,s3
void forward_pass(double* x)
{
	int i,j;
	double s;
	
	//perasma nevrwnwn 1ou epipedou
	for(i=0; i<H1; i++)
	{
		s=0;
		
		for(j=0; j<d; j++)
			s=s+w1[i][j]*x[j]; // athroisma x*w(eswterikou ginomenou)

		s1[i]=s;
		s=s+b1[i];
		z1[i]=f(s);	//eksodoi nevronwn prwtou epipedou (xrisi synartisi energopoihshs)		
		
	}
	
	
	//perasma nevronwn 2ou epipedou
	for(i=0; i<H2; i++)
	{
		s=0;
		
		for(j=0; j<H1; j++)
			s=s+w2[i][j]*z1[j]; // athroisma x*w(eswterikou ginomenou)

		s2[i]=s;
		s=s+b2[i];
		z2[i]=f(s);	//eksodoi nevronwn 2ou epipedou (xrisi synartisi energopoihshs)		
		
	}

	
	//perasma epipedou K(eksodou)
	for(i=0; i<K; i++)
	{
		s=0;
		
		for(j=0; j<H2; j++)
			s=s+w3[i][j]*z2[j]; // athroisma x*w(eswterikou ginomenou)

		s3[i]=s;
		s=s+b3[i];
		y[i]=f(s);	//eksodoi nevronwn 2ou epipedou (xrisi synartisi energopoihshs)		
		
	}
}


//backpropagation
void backprop(double* t)
{
	int i,j, sum;

	//sflama nevronwn eksodou
	//if(flag==1)
		for(i=0; i<K; i++)
			d3[i]=(y[i]-t[i])*y[i]*(1-y[i]);
	//f(flag==2)
		for(i=0; i<K; i++)
			d3[i]=(y[i]-t[i]);

	//sfalma 2ou epipedou
	for(i=0; i<H2; i++)
	{
		sum=0;
		for(j=0; j<K; j++)
			sum=sum+d3[j]*w3[j][i];
	
		d2[i]=sum*z2[i]*(1-z2[i]);
	}

	//sfalma 1ou epipedou
	for(i=0; i<H1; i++)
	{
		sum=0;
		for(j=0; j<H2; j++)
			sum=sum+d2[j]*w2[j][i];
	
		d1[i]=sum*z1[i]*(1-z1[i]);
	}
	
		
}

//ypologismos athoismatos paragwgon varwn
void sum_dw(double* x)
{
	int i,j;

	for(i=0; i<H1; i++)
		for(j=0; j<d; j++)
			dw1[i][j]=dw1[i][j]+(d1[i]*x[j]);

	for(i=0; i<H2; i++)
		for(j=0; j<H1; j++)
			dw2[i][j]=dw2[i][j]+(d2[i]*z1[j]);

	for(i=0; i<K; i++)
		for(j=0; j<H2; j++)
			dw3[i][j]=dw3[i][j]+(d3[i]*z2[j]);


}


//ypologismos athoismatos paragwgon polosewn
void sum_db()
{
	int i,j;

	for(i=0; i<H1; i++)
		db1[i]=db1[i]+d1[i];

	for(i=0; i<H2; i++)
		db2[i]=db2[i]+d2[i];

	for(i=0; i<K; i++)
		db3[i]=db3[i]+d3[i];


}




//gradient-descent
void gradient_descent()
{
	int i,j;

	//enimerwsi varwn+polosewn eksodou
	for(i=0;i<K;i++)
	{
		b3[i]=b3[i]-n*db3[i];

		for(j=0; j<H2; j++)
			w3[i][j]=w3[i][j]-n*dw3[i][j];

	}		

	//enimerwsi varwn+polosewn 2ou epipedou
	for(i=0;i<H2;i++)
	{
		b2[i]=b2[i]-n*db2[i];

		for(j=0; j<H1; j++)
			w2[i][j]=w2[i][j]-n*dw2[i][j];
	}

	//enimerwsi varwn+polosewn 1ou epipedou
	for(i=0;i<H1;i++)
	{
		b1[i]=b1[i]-n*db1[i];

		for(j=0; j<d; j++)
			w1[i][j]=w1[i][j]-n*dw1[i][j];
	}

}
//read from input files
void read_arxeio(){
  
    int i,j;
    double tempx1,tempx2;
    int sum=0;
    
    //protypa synolou ekpaideushs kai apothikeush tous se pinakes (prosoxi dimiourgountai mia fora)
    for(i=0;i<N;i++)
    {
	tempx1=GetRand();
	tempx2=GetRand();

	x_train[i][0]=tempx1;
	x_train[i][1]=tempx2;
	

	sum=pow(tempx1,2)+pow(tempx2,2);
	
	if(sum<=0.16)
	{
		
		 t_train[i][0]=0;
		 t_train[i][1]=0;
		 t_train[i][2]=1;
	}

	if(sum>0.16 && sum<=0.64)
	{
		
		 t_train[i][0]=0;
		 t_train[i][1]=1;
		 t_train[i][2]=0;
	}

	if(sum>0.64)
	{
		
		 t_train[i][0]=1;
		 t_train[i][1]=0;
		 t_train[i][2]=0;
	}
    }	

    //protypa synolou elegxou
    for(i=0;i<T;i++)
    {
	tempx1=GetRand();
	tempx2=GetRand();

	x_test[i][0]=tempx1;
	x_test[i][1]=tempx2;
	

	sum=pow(tempx1,2)+pow(tempx2,2);
	
	if(sum<=0.16)
	{
		
		 t_test[i][0]=0;
		 t_test[i][1]=0;
		 t_test[i][2]=1;
	}

	if(sum>0.16 && sum<=0.64)
	{
		
		 t_test[i][0]=0;
		 t_test[i][1]=1;
		 t_test[i][2]=0;
	}

	if(sum>0.64)
	{
		
		 t_test[i][0]=1;
		 t_test[i][1]=0;
		 t_test[i][2]=0;
	}
    }	


}



//sinartisi termatismou
int finish(){
    int i,j;
    double diafora,Ei;
    
    E=0;
    
    for(i=0;i<N;i++){
	for(j=0;j<K;j++){
	  
	      forward_pass(x_train[i]);
	      Ei=y[j]-t_train[i][j];
	      E=E+(((double)1/(double)2)*pow(Ei,2));
	}
      
      
    }
  //briskoume thn diafora metaksi prin sfalmatosk ai twrinou
  diafora=fabs(Epre-E);
  
  printf("\t\t  errors of the previous epoche %lf\n", Epre);
  printf("\t\t error of the current epoche %lf\n\n", E);
  Epre=E;
  
  //elenxoume
  if(diafora<e)
      return(1);
  else
      return(0);
  
}



//////------------------main--------

main()
{
  
  int i,j,epoxh=0;
	
    
  
 //tipwma eisodwn
 read_arxeio();
 
 
  
  //arxikopoihsh varwn-polosewn sto -1,1
  initialize_poloseis();
  initialize_varwn();
  
  //arxikopoihshs paragwgwn sto 0,0
  initialiaze_parag_polosewn();
  initialiaze_parag_varwn();
  
  //arxi epoxwn
  while(epoxh<MaxEpoxh)
  {
      for(i=0;i<N;i++)
      {
	
	  //euthi perasma
	  forward_pass(x_train[i]);
	  
	  //anapodo perasma
	  backprop(t_train[i]);
	  
	  //upologismos athroismatos paragogwn varwn-polosewn
	  sum_dw(x_train[i]);
	  sum_db();
	  
      }
    
      //gradient descent
      gradient_descent();
      
    
	    
	    
      //ta barh-polwseis epoxis
      //printf_varoi_polwseis();
      
	    
	    
    //kritirio termatismou
    if(finish()==1)
    {
	    printf("MLP is trained in:  %d \n",epoxh);
	    break;
    }
    //an oxi
    else
    {
	  epoxh++;
	  initialiaze_parag_varwn();
	  initialiaze_parag_polosewn();
    
    }
   
  
  }



	// vevaiotita 0.8
	//double vev=0.8;
	double sin1=0;
	double sin2=0;

	for(i=0; i<T; i++)
	{
		printf("%f  %f--->\n", x_test[i][0], x_test[i][1]);
		forward_pass(x_test[i]);
			
		for(j=0; j<K; j++)
		{
			printf("%f it should be %f\n", y[j], t_test[i][j]);
			if(y[j]>0.5)
			{
				sin1++;				
			}
			if(y[j]<0.5)
			{
				sin2++;
			}
			
		}
		
	}

	double sin=sin1+sin2;
	//printf("%lf\n", (sin/T)*100);
  	
	
}

