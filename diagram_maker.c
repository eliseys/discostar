#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "def.h"

int make_diagram(int dia)
{
  
  if (dia == 1)
    {
      FILE *diagram;
      diagram = fopen("DIAGRAM_JI_60", "a");

      int size_diagram = 180*20;
      int i, k, j;

      double PSI_pr = 0.0 * (M_PI/180.0);

      /* double * test_diagram = (double *) malloc(sizeof(double) * 180); */
      
      /* for (i = 0; i <180; i++) */
      /* 	{ */
      /* 	  if (i == 10 || i == 45 || i == 80 || i == 90 || i == 100 || i == 135 || i == 170) */
      /* 	    { */
      /* 	      test_diagram[i] = 1.0; */
      /* 	    } */
      /* 	  else */
      /* 	    { */
      /* 	      test_diagram[i] = 0.0; */
      /* 	    } */
      /* 	} */

      // i here is steps in PSI_pr 



      /* double * Ix_dd2; */

      /* printf("PSI_pr %f\n", PSI_pr); */
      
      /* Ix_dd2 = x_ray_direction_diagram(PSI_pr); */
	
      
      for (i = 0; i < 20; i++)
      	{
      	  PSI_pr = (double) i * (360.0/20.0) * (M_PI/180.0);

      	  double * Ix_dd2;


      	  Ix_dd2 = x_ray_direction_diagram(PSI_pr);
      	  //Ix_dd2 = test_diagram;


      	  /* for (j = 0; j < 180; j++) */
      	  /* 	{ */
      	  /* 	  printf("%f\n", Ix_dd2[j]); */
      	  /* 	} */
      
      	  fwrite(Ix_dd2, sizeof(double), 180, diagram);

      	  /* for (k = 0; k < 180; k++) */
      	  /* 	{ */
      	  /* 	  fprintf(diagram, "%.20f\t", Ix_dd2[k]); */

      	  /* 	} */
      
      	  //free(Ix_dd2);

      	  printf("PROGRESS %d\t %f\t %f\n", i, PSI_pr*(180.0/M_PI), Ix_dd2[1]);
      	}

      fclose(diagram);
  
      double test[size_diagram];

      FILE *data;
      data = fopen("DIAGRAM_JI_60", "r");
  
      fread(&test, sizeof(double), size_diagram, data);

      for (j = 0; j < 180*20; j++)
      	{
	  if (isnan(test[j]))
	    {
	      printf("BAD DIAGRAM");
	    }
	  
      	}
  
      fclose(data);

      
    }

  return 0;
}
