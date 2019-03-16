#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main()
{
  int diagram_size = 180*360;
  double diagram[diagram_size];

  FILE *file;
  file = fopen("DIAGRAM", "r");
  fread(&diagram, sizeof(double), diagram_size, file);
  fclose(file);


  int PSI_pr = 0;
  int PSI_pr_i;
  int i = 0;
  
  double *Ix_dd;

  Ix_dd = &diagram[0];

  double SUM1, SUM2;
  
  for(PSI_pr = 0; PSI_pr < 360; PSI_pr++)
    {
      int PSI_pr_i = (int) (360.0 * 180.0 * PSI_pr)/(360 * M_PI);

      SUM1 = 0.0;
      SUM2 = 0.0;
      for(i = 80; i <= 100; i++){

	if (i < 90){
	  SUM1 = SUM1 + Ix_dd[180*PSI_pr + i];
	}
	else if (i > 90){
	  SUM2 = SUM2 + Ix_dd[180*PSI_pr + i];
	}
	  
	//printf("%d\t%d\t%f\n", PSI_pr, i, Ix_dd[180*PSI_pr + i]);
	
      }
      printf("%d\t%f\t%f\n", PSI_pr, SUM1, SUM2);
      
    }



  return 0;
}
