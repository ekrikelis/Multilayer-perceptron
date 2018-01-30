#include <time.h>
#include <stdio.h>
#include <stdlib.h> 

//int init=0;

int init=0;
//dimiourgia tyxaiou arithmou
double GetRand()
{
 int i;
 while(init==0)
 {
    srand((unsigned)(time(0)));
    init=1;
 }
 
	  double rr = (   ((double)rand() / ((double)(RAND_MAX))*18+1.0));
	      
	  return(rr);
 
}





main()
{
		double rr;
	    int i;
	  //  initrand();
	    for(i=0;i<50;i++)
	    {
	       rr = (int)GetRand();
	      //double x=rr*M;
	      printf("%lf\n",rr);
	    }
}



/*

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main ()
{
  double iSecret, iGuess;


  srand ( time(NULL) );

  
  iSecret = rand() % 10 + 1;

  do {
    printf ("Guess the number (1 to 10): ");
    scanf ("%d",&iGuess);
    if (iSecret<iGuess) puts ("The secret number is lower");
    else if (iSecret>iGuess) puts ("The secret number is higher");
  } while (iSecret!=iGuess);

  puts ("Congratulations!");
  return 0;
}
*/
