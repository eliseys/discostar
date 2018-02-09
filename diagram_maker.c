#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "def.h"

int make_diagram(int dia)
{
  
  if (dia == 1)
    {
      FILE *diagram;
      diagram = fopen("DIAGRAM", "a");

      int size_diagram = 180*360;
      int i, k, j;

      double PSI_pr;
      
      for (i = 0; i < 360; i++)
      	{
      	  PSI_pr = (double) i * M_PI/180.0;

      	  double * Ix_dd2;

      	  Ix_dd2 = x_ray_direction_diagram(PSI_pr);

      	  /* for (j = 0; j < 180; j++) */
      	  /* 	{ */
      	  /* 	  printf("%f\n", Ix_dd2[j]); */
      	  /* 	} */
      
      	  fwrite(Ix_dd2, sizeof(double), 180, diagram);

      	  /* for (k = 0; k < 180; k++) */
      	  /* 	{ */
      	  /* 	  fprintf(diagram, "%.20f\t", Ix_dd2[k]); */

      	  /* 	} */
      
      	  free(Ix_dd2);

      	  printf("PROGRESS %d\n", i);
      	}

      fclose(diagram);
  
      double test[size_diagram];

      FILE *data;
      data = fopen("DIAGRAM", "r");
  
      fread(&test, sizeof(double), size_diagram, data);

      for (j = 0; j < 180*360; j++)
      	{
      	  printf("%f\t", test[j]);
      	}
  

      fclose(diagram);

      printf("Hello from dia!\n");
    }

  return 0;
}
