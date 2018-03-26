/* main file with example of how to run ECOS */
#include <stdio.h>
#include "ecos/ecos.h"
#include "ecos/data.h"

int main(void)
{
	/*char ver[7];*/
    idxint exitflag = ECOS_FATAL;
	pwork* mywork;

	/* set up data */	
	mywork = ECOS_setup(n, m, p, l, ncones, q, 0, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b);
 
    if( mywork != NULL )
    {
		/* solve */	
		exitflag = ECOS_solve(mywork);

		double var[n];
		printf("The solution is \n");
		for (int i = 0; i < n; i++)
        {
            var[i] = mywork->x[i];
            printf("%f", var[i]);
        }

    	/* clean up memory */
		ECOS_cleanup(mywork, 0);
    }
    
    /* explicitly truncate exit code */
	return (int)exitflag;
}
