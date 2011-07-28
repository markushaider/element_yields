#include <stdio.h>
#include <math.h>
#include <stdlib.h>


double  get_yields(double z, double m, double yields[][7]);

const double ytable_z[4]={0,0.001,0.004,0.02};
const double ytable_m[7]={13,15,18,20,25,30,40};

/* SN II Yields from Nomoto 2006 */
/* Nucleosynthesis Yields of Core-Collapse Supernovae and */
/* Hypernovae, and Galactic Chemical Evolution */
double yields_fe[4][7]={
  {7.00E-2, 7.00E-2, 7.00E-2, 7.00E-2, 7.00E-2, 7.00E-2, 7.02E-2},
  {7.26E-2, 7.08E-2, 7.11E-2, 7.09E-2, 7.11E-2, 7.12E-2, 7.15E-2},
  {7.26E-2, 7.30E-2, 8.72E-2, 7.40E-2, 7.47E-2, 7.46E-2, 7.47E-2},
  {8.32E-2, 8.51E-2, 8.72E-2, 8.87E-2, 9.01E-2, 9.18E-2, 8.08E-2}
};

double yields_ox[4][7]={
  {4.50E-1, 7.73E-1, 1.38   , 2.11   , 2.79, 4.81, 8.38},
  {5.04E-1, 2.94E-1, 4.22E-1, 2.18   , 3.82, 5.33, 8.37},
  {3.85E-1, 2.92E-1, 5.21E-1, 9.94E-1, 2.20, 4.79, 7.96},
  {2.18E-1, 1.62E-1, 7.70E-1, 1.05   , 2.35, 3.22, 7.33}
};

double yields_si[4][7]={
  {8.04E-2, 7.32E-2, 1.16E-1, 9.94E-2, 3.51E-1, 2.48E-1, 1.02   },
  {8.99E-2, 4.29E-2, 1.53E-1, 1.28E-1, 1.20E-1, 1.65E-1, 8.81E-1},
  {6.11E-2, 1.03E-1, 9.41E-2, 1.24E-1, 1.19E-1, 3.95E-1, 5.23E-1},
  {7.48E-2, 8.38E-2, 1.01E-1, 6.32E-2, 1.28E-1, 2.40E-1, 2.41E-1}
};

/* sum of all metals from Nomoto */
double yields_all[4][7]={
  {0.82, 1.54, 2.50, 3.63, 4.41, 6.71, 11.19},
  {0.98, 0.79, 1.13, 3.49, 5.73, 7.57, 10.88},
  {0.85, 0.83, 1.45, 1.78, 3.70, 6.94, 11.71},
  {0.67, 0.60, 1.54, 2.11, 4.30, 5.38, 11.36}
};

/* W7 */
double snIa_yields_fe =0.626;
double snIa_yields_ox =0.143;
double snIa_yields_si =0.154;

/* WDD2 */
/* double snIa_yields_fe =0.713 */
/* double snIa_yields_ox =0.0658 */
/* double snIa_yields_si =0.206 */

int main(int argc, char** argv) {

  char *lpath, *spath;
  lpath = argv[1];
  spath = argv[2];

  /* Berechne Tabelle für die IMF */
  const double alpha[3] = {-0.3,-1.3,-1.7};
  /* Normierung für den Bereich m = 0.08 - 100.0 */
  double norm = pow(1.0/0.08,alpha[0])/(alpha[0]+1)*(pow(0.5,alpha[0]+1)-pow(0.08,alpha[0]+1))
    + pow(0.5/0.08,alpha[0])*pow(1/0.5,alpha[1])/(alpha[1]+1)*(pow(1.0,alpha[1]+1)-pow(0.5,alpha[1]+1))
    + pow(0.5/0.08,alpha[0])*pow(1/0.5,alpha[1])/(alpha[2]+1)*(pow(100.0,alpha[2]+1)-pow(1.0,alpha[2]+1));
  double imf_kroupa[1000];
  /* Imf ausrechnen */
  int i;
  for(i=0;i<1000;i++) {
    if((double)i/10.0<0.5)
      imf_kroupa[i]=pow((double)i/(10.0*0.08),alpha[0])/norm;
    else if(0.5<=(double)i/10.0 && (double)i/10.0<1.0)
      imf_kroupa[i]=pow(0.5/0.08,alpha[0])*pow((double)i/(10.0*0.5),alpha[1])/norm;
    else if((double)i/10.0>=1.0)    
      imf_kroupa[i]=pow(0.5/0.08,alpha[0])*pow(1.0/0.5,alpha[1])*pow((double)i/10.0,alpha[2])/norm;
  }
  /* Achtung, imf_kroupa[0] ist inf, da Hochzahl negativ */
  /* Debug Output */
  /* i=0; */
  /* for(i;i<500;i++) { */
  /*   printf("%g %g\n",(double)i/10.0,imf_kroupa[i]);  */
  /* } */

  /* Einlesen der galaxy files */
  /* Number of galaxy files, equals number of lines in file "galaxies" */
  int nfiles=1286;
  /* Number of galaxies (=lines) in each galaxy file */
  int ngalxs=36268;
  /* Array für Galaxy-Dateinamen und Zeit */
  char **fnames = (char **)malloc(sizeof(char*)*nfiles);
  double *time  = (double *)malloc(sizeof(double)*nfiles);
  double *line  = (double *)malloc(sizeof(double)*48); //(48 fields)
 
  
  /* Öffnen von "galaxies" welches Dateinamen und Zeit beinhaltet */
  char filename[160];
  sprintf(filename,"%s/galaxies",lpath);
  FILE *gfile = fopen(filename,"r");
  if(gfile == NULL)
  {
    printf("%s\n",filename);
    printf("Cannot open galaxies file\n");
    exit(-1);
  }

  /* Einlesen der Dateinamen und der Zeit */  
  i=-1;
  do {
    i++;
    fnames[i] = (char*)malloc(sizeof(char)*15);
  }
  while(fscanf(gfile,"%s %*f %*f %*f %lf %*i\n",fnames[i],&time[i])==2); 
  fclose(gfile);
  for (i =0; i<nfiles; i++) 
  printf("name %s, time %lf\n",fnames[i],time[i]);

  char fname[100];
  //double *stf_tot = (double *)calloc(nfiles*,sizeof(double));
  //double stf[nfiles][ngalxs];

  double **stf = (double **)malloc(nfiles * sizeof(double *));
  for(i = 0; i < nfiles; i++)
    stf[i] = (double *)calloc(ngalxs,sizeof(double));
  double **rdisk = (double **)malloc(nfiles * sizeof(double *));
  for(i = 0; i < nfiles; i++)
    rdisk[i] = (double *)calloc(ngalxs,sizeof(double));
  double **metal = (double **)malloc(nfiles * sizeof(double *));
  for(i = 0; i < nfiles; i++)
    metal[i] = (double *)calloc(ngalxs,sizeof(double));
  double **init_metal = (double **)malloc(nfiles * sizeof(double *));
  for(i = 0; i < nfiles; i++)
   init_metal[i] = (double *)calloc(ngalxs,sizeof(double));


  /* Schleife über die galaxy files */
  int gfn,ngal,k;
  /* Zeitschritt 0 außerhalb der Schleife */
  sprintf(fname,"%s/%s",lpath,fnames[0]);
  gfile = fopen(fname,"r");
  if(gfile == NULL) {
    printf("Cannot open file %s\n",fnames[0]);
    exit(-1);
  }
  for (ngal=0; ngal<ngalxs;ngal++) {
    for (k = 0;k<48; k++) {
      fscanf(gfile,"%lf",&line[k]);
    }
    rdisk[0][ngal]=line[0];
    if (line[0]>0&&line[2]>1E-3) {
      stf[0][ngal]=line[33]+line[34]+line[35]+line[36]+line[37];
      metal[0][ngal]=(line[38]+line[39]+line[40]+line[41]+line[42])/5*0.02;
      init_metal[0][ngal]=line[2]*1.0E10*line[19];
    }
  }
  fclose(gfile);
  for (gfn = 1; gfn<nfiles; gfn++)	
  {
    printf("fname %s\n",fnames[gfn]);
    sprintf(fname,"%s/%s",lpath,fnames[gfn]);
    gfile = fopen(fname,"r");
    if(gfile == NULL) {
      printf("Cannot open file %s\n",fnames[gfn]);
      exit(-1);
    }
    /* Einlesen der galaxy Dateien */
    for( ngal = 0; ngal<ngalxs; ngal++) {
      for (k = 0;k<48; k++) {
 	fscanf(gfile,"%lf",&line[k]);
      }
      fscanf(gfile,"\n");
      rdisk[gfn][ngal]=line[0];
      if (line[0]>0&&line[2]>1E-3) {
	stf[gfn][ngal]=line[33]+line[34]+line[35]+line[36]+line[37];
	metal[gfn][ngal]=(line[38]+line[39]+line[40]+line[41]+line[42])/5.0*0.02;
	if (rdisk[gfn-1][ngal]==0) {
	  /* nach Anders, Grevesse 1989 */
	  /* in solar masses, neu: falsch, da man noch mit 0.02 multiplizieren müsste */
	  init_metal[gfn][ngal]=line[2]*1.0E10*line[19];
	}
      }
    }
    fclose(gfile);
  }
  free(line);
  for(i = 0; i < nfiles; i++)
    free(rdisk[i]);
  free(rdisk);

  double **snrIa = (double **)malloc(nfiles * sizeof(double *));
  for(i = 0; i < nfiles; i++)
    snrIa[i] = (double *)calloc(ngalxs,sizeof(double));
  double **snrII = (double **)malloc(nfiles * sizeof(double *));
  for(i = 0; i < nfiles; i++)
    snrII[i] = (double *)calloc(ngalxs,sizeof(double));
  double **snII_fe = (double **)malloc(nfiles * sizeof(double *));
  for(i = 0; i < nfiles; i++)
    snII_fe[i] = (double *)calloc(ngalxs,sizeof(double));
  double **snII_ox = (double **)malloc(nfiles * sizeof(double *)); 
  for(i = 0; i < nfiles; i++)
    snII_ox[i] = (double *)calloc(ngalxs,sizeof(double));
  double **snII_si = (double **)malloc(nfiles * sizeof(double *));
  for(i = 0; i < nfiles; i++)
    snII_si[i] = (double *)calloc(ngalxs,sizeof(double));
  double **snIa_fe = (double **)malloc(nfiles * sizeof(double *));
  for(i = 0; i < nfiles; i++)
    snIa_fe[i] = (double *)calloc(ngalxs,sizeof(double));
  double **snIa_ox = (double **)malloc(nfiles * sizeof(double *));
  for(i = 0; i < nfiles; i++)
    snIa_ox[i] = (double *)calloc(ngalxs,sizeof(double));
  double **snIa_si = (double **)malloc(nfiles * sizeof(double *));
  for(i = 0; i < nfiles; i++)
    snIa_si[i] = (double *)calloc(ngalxs,sizeof(double));
  double **snII_metal = (double **)malloc(nfiles * sizeof(double *));
  for(i = 0; i < nfiles; i++)
   snII_metal[i] = (double *)calloc(ngalxs,sizeof(double));

  /* Füllen der SN Arrays */

  /* Berechne den Anteil der TypII Supernova */
  double integrated_imf=0;
  for(i=80;i<1000;i++) {
    integrated_imf+=imf_kroupa[i]*0.1;
  }
  int m;
  for(gfn=0;gfn<nfiles;gfn++) {
    for(ngal=0;ngal<ngalxs;ngal++) {
      snrII[gfn][ngal]=stf[gfn][ngal]*integrated_imf; /* m_sol pro Jahr */
      if(gfn+100<nfiles)
	snrIa[gfn+100][ngal]=stf[gfn][ngal]*0.022*0.02;	/* events pro jahr */
      for(m=80;m<400;m+=4) {
	/* synthesized solar masses */
	/* check again if necessary to divide by m */
	snII_fe[gfn][ngal]+=stf[gfn][ngal]*get_yields(metal[gfn][ngal],(double)m/10.0,yields_fe)*imf_kroupa[m]*0.4/m*10;
	snII_ox[gfn][ngal]+=stf[gfn][ngal]*get_yields(metal[gfn][ngal],(double)m/10.0,yields_ox)*imf_kroupa[m]*0.4/m*10;
	snII_si[gfn][ngal]+=stf[gfn][ngal]*get_yields(metal[gfn][ngal],(double)m/10.0,yields_si)*imf_kroupa[m]*0.4/m*10;
	snII_metal[gfn][ngal]+=stf[gfn][ngal]*get_yields(metal[gfn][ngal],(double)m/10.0,yields_all)*imf_kroupa[m]*0.4/m*10;
      }
    }
  }

  for(gfn=0;gfn<nfiles;gfn++) {
    for(ngal=0;ngal<ngalxs;ngal++) {
      /* synthesized solar masses */
      snIa_fe[gfn][ngal]=snrIa[gfn][ngal]*snIa_yields_fe;
      snIa_ox[gfn][ngal]=snrIa[gfn][ngal]*snIa_yields_ox;
      snIa_si[gfn][ngal]=snrIa[gfn][ngal]*snIa_yields_si;
    }
  }

  /* Ausschreiben */
  for (gfn = 0; gfn<nfiles; gfn++) {
    printf("fname %s\n",fnames[gfn]);
    sprintf(fname,"%s/%s",spath,fnames[gfn]);
    gfile = fopen(fname,"w");
    if(gfile == NULL) {
      printf("Cannot open file %s\n",fnames[gfn]);
      exit(-1);
    }
    for(ngal = 0; ngal<ngalxs; ngal++) {
      fprintf(gfile,"%.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e\n",
	      stf[gfn][ngal],metal[gfn][ngal],snrIa[gfn][ngal],snrII[gfn][ngal],snII_fe[gfn][ngal],
	      snII_ox[gfn][ngal],snII_si[gfn][ngal],snIa_fe[gfn][ngal],snIa_ox[gfn][ngal],snIa_si[gfn][ngal],
	      init_metal[gfn][ngal], snII_metal[gfn][ngal]);
    }
    fclose(gfile);
  }

  for(i = 0; i < nfiles; i++)
    free(stf[i]);
  free(stf);
  for(i = 0; i < nfiles; i++)
    free(metal[i]);
  free(metal);
  for(i = 0; i < nfiles; i++)
    free(snrIa[i]);
  free(snrIa);
  for(i = 0; i < nfiles; i++)
    free(snrII[i]);
  free(snrII);
  for(i = 0; i < nfiles; i++)
    free(snIa_fe[i]);
  free(snIa_fe);
  for(i = 0; i < nfiles; i++)
    free(snIa_ox[i]);
  free(snIa_ox);
  for(i = 0; i < nfiles; i++)
    free(snIa_si[i]);
  free(snIa_si);
  for(i = 0; i < nfiles; i++)
    free(snII_fe[i]);
  free(snII_fe);
  for(i = 0; i < nfiles; i++)
    free(snII_ox[i]);
  free(snII_ox);
  for(i = 0; i < nfiles; i++)
    free(snII_si[i]);
  free(snII_si);
  for(i = 0; i < nfiles; i++)
    free(snII_metal[i]);
  free(snII_metal);
  for(i = 0; i < nfiles; i++)
    free(init_metal[i]);
  free(init_metal);
  free(time);
  for(i = 0; i < nfiles; i++)
    free(fnames[i]);
  free(fnames);

  return 0;
}


double get_yields(double z,double m, double yields[][7]) {
  int i, j;
  for(i=0;i<4-1;i++){
    if(z<ytable_z[i])
      break;
  }
  for(j=0;j<7-1;j++) {
    if(m<ytable_m[j])
      break;
  }
  if(i==0)
    i=1;
  if(j==0)
    j=1;
  double p1,p2;
  p1=yields[i-1][j-1]+(yields[i][j-1]-yields[i-1][j-1])*(z-ytable_z[i-1])/(ytable_z[i]-ytable_z[i-1]);
  p2=yields[i-1][j]  +(yields[i][j]  -yields[i-1][j])  *(z-ytable_z[i-1])/(ytable_z[i]-ytable_z[i-1]);
  return p1+(p2-p1)*(m-ytable_m[j-1])/(ytable_m[j]-ytable_m[j-1]);
}
