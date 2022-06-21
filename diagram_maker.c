#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "def.h"

int main(){
  
  FILE *diagram;
  diagram = fopen("DIAGRAM_50", "a");

  int psi_step_num = 360; /* number of steps along precession axis*/
  int theta_step_num = 180; /* number of theta steps */
  
  int size_diagram = theta_step_num*psi_step_num; /* theta * phi */
  int i, j;
  double PSI_pr;
  
  for (i = 0; i < psi_step_num; i++)
    {
      PSI_pr = (double) i * M_PI/180.0;
      double * Ix_dd;

      Ix_dd = x_ray_direction_diagram(PSI_pr);
      
      fwrite(Ix_dd, sizeof(double), theta_step_num, diagram);
      
      free(Ix_dd);

      printf("PROGRESS %d\n", i);
    }

  fclose(diagram);
  


  double test[size_diagram];

  FILE *data;
  data = fopen("DIAGRAM_50", "r");
  
  fread(&test, sizeof(double), size_diagram, data);

  for (j = 0; j < theta_step_num; j++)
    {
      printf("%f\t", test[theta_step_num + j]);
    }
  
  
  //free(Ix_dd);

  fclose(diagram);

  return 0;
}
