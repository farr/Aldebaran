
/****************************************************/
/*	param.c	*/
/****************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
//#include <iostream> //TSR: testing dynamic vectors
//#include <vector> //TSR: testing dynamic vectors

#include "parametri.h"
#include "colors.h"
#include "tracks.h"
#include "sismo.h"
#include "pdf.h"
//#include "star.h"
#define	HR_RES	50 

/********************** estima from seismology ************************/
/* From an input table with Teff, [M/H], (numax and/or Dnu and/or and/or DP and/or luminosity),
 computes stellar properties (STEP 1): age, mass, radius, logg, mean density;
 adding apparent magnitudes computes distances and extinctions simultaneously */

int metododireto;
float use_teff, use_dnu, use_numax, use_lg, use_pspacing, use_lum,
  max_smagobs, delta_teff;
void alamb_av_teff (float, float *);
void read_param_conf ();

void	estima_from_seismology (isoc)
     ISOCRONA	*isoc ;
{

  //#define VERBOSE
#define NPROB_MU 3201 //TSR 3000 -> 3201
#define NPROB_JPDF 200000 //TSR 80000 -> 95000
#define NPROB_AV 901 //TSR 151 -> 351 -> 901
#define NSIGMA 3.0
#define ANNUL -99.9
  //#define SCRIVE_ISOC

  int i, j, ii, iii, nage, n_int,
    imu, iav, evstate,
    ic, nic, ijpdf;   
  float age, logl, logte, logg, logr,
    teff0, feh, fehv, offset_feh, delta_feh, sfeh, 
    feh_min, feh_max,
    mv0, r0, lrho, deltanu, numax0,
    mass_ini, mass_act, zeta,
    lteff, slteff,
    delta_mass_ini, delta_mass_act, delta_logl, delta_logte,
    weight, nstelle, fact_age, fact_feh,
    agemin, agemax,
    lage, lage_ini, lage_fin, delta_lage;
  float pspacing, delta_pspacing, delta_dnu;  // only for seismic tracks: MESA
  double lum, elum, log_lum, slog_lum;
  float lg, slg, pspac, spspac; 
  double dnu, sdnu, numax, snumax, mh, smh, teff, steff,
    gc_l, gc_b, avx, savx;  
  int fit_av, corr_mag, bright_star;
  double logrhodir, elogrhodir, massdir, emassdir, logmassdir, elogmassdir, 
    raiodir, eraiodir, lograiodir, elograiodir, loggdir, eloggdir;
  float emag;
  MAGNITU  deltamag;
  float mu0_av_prob[NPROB_JPDF][3],
    jpdf_mu0_av[NPROB_AV][NPROB_MU],
    mu0_min, mu0_max, mu0_mode, av_mode;
  long double volume; 
  MAGNITU  magni, magobs, smagobs, mu[NPROB_MU];
  float loggrav(), imf_n(), gaussiana();
  void	fa_isocrona_completa() ;
  char  header[2000], inpend[1000], inputfile[80], outputfile[80], temp[65];
  char  name[20], nomimag[40];
  clock_t start, end;
  time_t timer;
  struct tm ltime;

  PDF pdf_age, pdf_m, pdf_g, pdf_r, pdf_rho,
    pdf_mag[TAB_NCOL], pdf_magobs[TAB_NCOL], pdf_mu[TAB_NCOL],
    pdf_mu0_lamb[TAB_NCOL], pdf_av, pdf_mu0,
    pdf_loglum;
  PDF2 jpdf_age_m, jpdf_age_g, jpdf_age_r,
    jpdf_m_g, jpdf_m_r, jpdf_r_g;
  //Input PDFs
  PDF pdf_logte, pdf_fehv,
    pdf_deltanu, pdf_pspacing, pdf_numax0;
  
  //#define ANDREAOUTPUT
#ifdef ANDREAOUTPUT
  int im, ir;
  FILE *fp_out_tmp, *fppp_out;
  float prob[108][851][211],
    probtot;
#endif
 
  printf("\nFunction 32: Estimating distances/extinctions from (Dnu,nu_max,T_eff,[M/H]).");

#ifdef VERBOSE
  printf("\nGrid resolution = %d", HR_RES);
#endif

  /* Timer */
  char buffer[25];
  timer = time (0);
  localtime_r (& timer, & ltime);
  strftime (buffer, 25, "%A, %d/%m/%Y", & ltime);
  printf ("\n%s\n", buffer);
  
  start = clock();

  /* Defining photometric system for each sample -- NEED TO BE IMPROVED */
  char *fphotsys,
    photsys[4][11] = {"apokasc", "apogee", "corotsismo", "dat"} ;
  int ip; //apokasc = 0, apogee = 1, corotsismo = 2, others = 3
  int ext_coef_table ;
  for (ip=0;ip<4;ip++) 
    {
      fphotsys = strstr(file_mag, photsys[ip]);
      if (fphotsys != NULL)
	break;
    }
#ifdef VERBOSE
  printf("\n%s => Photometric system = %d -> %s\n", file_mag, ip, photsys[ip]);
#endif
  if (ip == 0) //APOKASC
    printf ("\nPhotometric system: combination of AB and VEGA. No correction need to be done.\n");
  else if (ip == 1 || ip == 2) //APOGEE_NEW e COROTSISMO
    printf ("\nPhotometric system: VEGA. Applying correction on ugriz apparent magnitudes...\n") ;
  else printf ("\nCheck if the input photometry is in the same photometric system as the models.");
  
  /* Setting extinction coeficient table */
  ext_coef_table = (ip == 0) ? 1 : 0; //only when using tab_mag_apokasc.dat
#ifdef VERBOSE
  ext_coef_table == 1 ?
    printf("\nUsing extinction coeficient table.\n") :
    printf("\nUsing Solar extinction coeficient table.\n") ;
#endif  

  //Reading PARAM configuration
  read_param_conf ();
  
  /* Creating PDFs directory */   
  if (print_pdfs==1)
    {
      strcpycat (temp, PDFSDIR, "PDFdist/");
      makedir (temp);
    }
  if (print_jpdfs==1)
    {
      strcpycat (temp, PDFSDIR, "JPDF_mag/");
      makedir (temp);
      strcpycat (temp, PDFSDIR, "JPDF_mu0_av/");
      makedir (temp);
    }
  if (ext_coef_table == 1)
    {
      strcpycat (temp, PDFSDIR, "AlambAv/");
      if (print_pdfs==1) makedir (temp) ;
    }
  else
    {
      struct stat st;
      char rmdir[30];
      strcpycat (temp, PDFSDIR, "AlambAv/");
      if(stat(temp, &st) == 0)
	{
	  strcpycat (rmdir, "rm -rf ", temp);
	  system (rmdir);
	}
    }
   
  //  vector<vector<double> > a;//define vector of vectors (matrix) TSR: trying to implement dynamic vectors  

  printf ("\n\t Inputfile = ");
  scanf ("%s", inputfile); 
  fp_inp = fopen(inputfile, "r");
  for (i=0; i<1; i++) //jumps header lines
    fgets(header, 2000, fp_inp) ;
  
  printf ("\n\t Median outputfile  = ");
  scanf ("%s", outputfile);
  fp_me_out = fopen(outputfile, "w");
  
  printf ("\n\t Mode outputfile  = ");
  scanf ("%s", outputfile);
  fp_mo_out = fopen(outputfile, "w");
  
  if (print_input_pdfs==1)
    {
      char me_inpoutputfile[80],  mo_inpoutputfile[80];
      printf ("\n\t Median input_PDF outputfile  = ");
      scanf ("%s", me_inpoutputfile);
      fp_inppdfs_me = fopen(me_inpoutputfile, "w");
      printf ("\n\t Mode input_PDF outputfile  = ");
      scanf ("%s", mo_inpoutputfile);
      fp_inppdfs_mo = fopen(mo_inpoutputfile, "w");
    }

  //Median output
  print_outputheader(fp_me_out, "Id age mass logg rad logrho");
  for (ic=0; ic<n_fil+1; ic++)
    print_outputheader(fp_me_out, nomi_mag[ic]);  
  print_outputheader(fp_me_out,"Av mu0 dist nfil fils");
  fprintf (fp_me_out, "\n");
  //Mode output
  print_outputheader(fp_mo_out, "Id age mass logg rad logrho");
  for (ic=0; ic<n_fil+1; ic++)
    print_outputheader(fp_mo_out, nomi_mag[ic]);  
  print_outputheader(fp_mo_out,"Av2d Av mu02d mu0 dist2d dist nfil fils");
  fprintf (fp_mo_out, "\n");
  //Input_pdfs output
  if (print_input_pdfs==1)
    {
      print_outputheader(fp_inppdfs_me, "Id Teff MH"); print_outputheader(fp_inppdfs_mo, "Id Teff MH");
      if (use_dnu>0.0)
	{
	  print_outputheader(fp_inppdfs_me, "Dnu"); print_outputheader(fp_inppdfs_mo, "Dnu");
	}
      if (kind_tracks==3 && use_pspacing>0.0)
	{
	  print_outputheader(fp_inppdfs_me, "DPi1"); print_outputheader(fp_inppdfs_mo, "DPi1");
	}
      if (use_numax>0.0)
	{
	  print_outputheader(fp_inppdfs_me, "numax"); print_outputheader(fp_inppdfs_mo, "numax");
	}
      if (use_lum>0.0)
        {
          print_outputheader(fp_inppdfs_me, "lum"); print_outputheader(fp_inppdfs_mo, "lum");
        }
      fprintf (fp_inppdfs_me, "\n"); fprintf (fp_inppdfs_mo, "\n");
    }

  offset_feh = 0.0 ;
  printf ("\n\t Min, maximum searched age(yr) = ");
  scanf ("%f %f", &agemin, &agemax);

  /* Input table:
     #ID Dnu eDnu numax enumax pspac epspacing evstate
     feh efeh teff eteff logg elogg l b
     mag1 emag1 mag2 emag2 mag3 emag3 ... magN emagN
     Av eAv
     additional columns that it will be not read by param 
  */
  while (fscanf (fp_inp,  
		 "%s %lf %lf %lf %lf %f %f  %d %lf %lf %lf %lf %f %f  %lf %lf", 
		 name, &dnu, &sdnu, &numax, &snumax, &pspac, &spspac,
		 &evstate, &mh, &smh, &teff, &steff, &lg, &slg,
  		 &gc_l, &gc_b) != EOF )
    {
      for (ic=1; ic<n_fil+1; ic++)
	fscanf(fp_inp, "%f %f", &magobs.fil[ic], &smagobs.fil[ic]) ;     
      fscanf(fp_inp, "%lf %lf", &avx, &savx); // READING Av HERE!
      if (use_lg > 0.0) printf ("\n Using logg=%.3f+-%.3f", lg, slg);
      if (use_dnu > 0.0) printf ("\n Using Dnu=%.3f+-%.3f muHz", dnu, sdnu);
      if (use_numax > 0.0) printf ("\n Using nu_max=%.3f+-%.3f muHz", numax, snumax);
      if (use_pspacing > 0.0) printf ("\n Using DP=%.3f+-%.3f muHz", pspac, spspac);
      if  (use_lum > 0.0)
	{
	  fscanf(fp_inp, "%lf %lf", &lum, &elum); // READING LUMINOSITY HERE!
	  log_lum = log10(lum);
	  slog_lum = elum / (lum * log(10.0)); ///log_lum-log10(lum-elum);
	  printf ("\n Using luminosity: L/Ls=%.3f+-%.3f, log(L/Ls)=%.3f+-%.3f", lum, elum, log_lum, slog_lum);
	}
      printf ("\n");
      fgets (inpend, 1000, fp_inp);
      
      /* Correcting small uncertainty in magnitudes */
      for (ic=1; ic<n_fil+1; ic++)
	if (smagobs.fil[ic] > -90.0 && smagobs.fil[ic] <= 0.01) smagobs.fil[ic]=0.011 ;
      
      /* Setting maximum uncertainty in magnitudes */ //For BRIGHT STARS use: 0.1 -> 0.9
      for (ic=1; ic<n_fil+1; ic++)
	if (smagobs.fil[ic] >= max_smagobs)
	  {
	    magobs.fil[ic] = ANNUL; smagobs.fil[ic] = ANNUL ;
	  }
      
      /* Setting null magnitude values and
	 correcting zero point ABmag -> Vegamag systems*/
      float u_cor = -0.925, g_cor = 0.107,
	r_cor = -0.142, i_cor = -0.355, z_cor = -0.518;
      
      if (ip == 0) /* tab_mag_apokasc.dat */
	{
	  magobs.fil[1] = ANNUL; smagobs.fil[1] = ANNUL; //Kep
	  magobs.fil[6] = ANNUL; smagobs.fil[6] = ANNUL; //d51
	  magobs.fil[12] = ANNUL; smagobs.fil[12] = ANNUL ; //W3
	  magobs.fil[13] = ANNUL; smagobs.fil[13] = ANNUL; //W4
	}
      else if (ip == 1) /* tab_mag_apogee_new.dat */
	{
	  magobs.fil[1] = ANNUL; smagobs.fil[1] = ANNUL ; //uS 
	  
	  //WITHOUT KIC and i(APASS/SDSS) - TSR
	  
	  //magobs.fil[2] = ANNUL; smagobs.fil[2] = ANNUL ; //gS
	  //magobs.fil[3] = ANNUL; smagobs.fil[3] = ANNUL ; //rS
	  magobs.fil[4] = ANNUL; smagobs.fil[4] = ANNUL ; //iS
	  //magobs.fil[5] = ANNUL; smagobs.fil[5] = ANNUL ; //zSx
	 
	  magobs.fil[6] = ANNUL; smagobs.fil[6] = ANNUL ; //gK
	  magobs.fil[7] = ANNUL; smagobs.fil[7] = ANNUL ; //rK
	  magobs.fil[8] = ANNUL; smagobs.fil[8] = ANNUL ; //iK
	  magobs.fil[9] = ANNUL; smagobs.fil[9] = ANNUL ; //zK
	  
	  //magobs.fil[10] = ANNUL; smagobs.fil[10] = ANNUL ; //gA
	  //magobs.fil[11] = ANNUL; smagobs.fil[11] = ANNUL ; //rA
	  magobs.fil[12] = ANNUL; smagobs.fil[12] = ANNUL ; //iA
	  //magobs.fil[13] = ANNUL; smagobs.fil[13] = ANNUL ; //BA
	  //magobs.fil[14] = ANNUL; smagobs.fil[14] = ANNUL ; //VA
	  
	  if (magobs.fil[1] > -90.0) magobs.fil[1] += u_cor; //uS
	  if (magobs.fil[2] > -90.0) magobs.fil[2] += g_cor; //gS
	  if (magobs.fil[3] > -90.0) magobs.fil[3] += r_cor; //rS
	  if (magobs.fil[4] > -90.0) magobs.fil[4] += i_cor; //iH    
	  if (magobs.fil[5] > -90.0) magobs.fil[5] += z_cor; //zH    
	  
	  if (magobs.fil[6] > -90.0) magobs.fil[6] += g_cor; //gK
	  if (magobs.fil[7] > -90.0) magobs.fil[7] += r_cor; //rK
	  if (magobs.fil[8] > -90.0) magobs.fil[8] += i_cor; //iK
	  if (magobs.fil[9] > -90.0) magobs.fil[9] += z_cor; //zK
	  
	  if (magobs.fil[10] > -90.0) magobs.fil[10] += g_cor; //gA
	  if (magobs.fil[11] > -90.0) magobs.fil[11] += r_cor; //rA
	  if (magobs.fil[12] > -90.0) magobs.fil[12] += i_cor; //iA 
	}
      else if (ip == 2) /* tab_mag_corotsismo.dat */
	{
	  magobs.fil[1] = ANNUL; smagobs.fil[1] = ANNUL ; //u_mag_sdss 
	  
	  //WITHOUT APASS - TSR
	  /*   
	       magobs.fil[7] = ANNUL; smagobs.fil[7] = ANNUL; //gA
	       magobs.fil[8] = ANNUL; smagobs.fil[8] = ANNUL; //rA
	       magobs.fil[9] = ANNUL; smagobs.fil[9] = ANNUL; //iA
	       magobs.fil[10] = ANNUL; smagobs.fil[10] = ANNUL; //BA
	       magobs.fil[11] = ANNUL; smagobs.fil[11] = ANNUL; //VA
	  */
	  
	  if (magobs.fil[1] > -90.0) magobs.fil[1] += u_cor; //uH
	  if (magobs.fil[2] > -90.0) magobs.fil[2] += g_cor; //gH
	  if (magobs.fil[3] > -90.0) magobs.fil[3] += r_cor; //rH
	  if (magobs.fil[4] > -90.0) magobs.fil[4] += i_cor; //iH      
	  if (magobs.fil[7] > -90.0) magobs.fil[7] += g_cor; //gA
	  if (magobs.fil[8] > -90.0) magobs.fil[8] += r_cor; //rA
	  if (magobs.fil[9] > -90.0) magobs.fil[9] += i_cor; //iA     
	  
	  /* Making sure that APASS is compatible with SLOAN-GUNN. CoRoGEE */
	  for (ic=7; ic<12; ic++)
	    {
	      if (magobs.fil[ic] <= -90.0 || magobs.fil[ic-5] <= -90.0) continue;
	      // printf ("\nAPASS %s %f %f SLOAN-GUNN %s %f %f", nomi_mag[ic], magobs.fil[ic], smagobs.fil[ic],
	      //	  nomi_mag[ic-5], magobs.fil[ic-5], smagobs.fil[ic-5]);
	      if ((magobs.fil[ic] + 3.0*smagobs.fil[ic] <= magobs.fil[ic-5]-3.0*smagobs.fil[ic-5]) ||
		  (magobs.fil[ic] - 3.0*smagobs.fil[ic] >= magobs.fil[ic-5]+3.0*smagobs.fil[ic-5]))
		{
		  printf ("\nAPASS %s (%f +- %f) is not compatible with SLOAN-GUNN %s  (%f +- %f)",
			  nomi_mag[ic], magobs.fil[ic], smagobs.fil[ic],
			  nomi_mag[ic-5], magobs.fil[ic-5], smagobs.fil[ic-5]);
		  magobs.fil[ic] = ANNUL; smagobs.fil[ic] = ANNUL;
		}
	    }
	}
      // printf ("\n");
 
      /* Printing input table */      
      printf("\n%s %f %f %f %f %f %f  %d %f %f %f %f %f %f   %f %f",
      	     name, dnu, sdnu, numax, snumax, pspac, spspac,
	     evstate, mh, smh, teff, steff, lg, slg, gc_l, gc_b) ;
      for (ic=1; ic<n_fil+1; ic++)	     
	printf(" %f %f ", magobs.fil[ic], smagobs.fil[ic]);
      printf ("%f %f", avx, savx); 
      if (use_lum > 0) printf ("%f %f", lum, elum);
      printf (" %s", inpend); 
      
      fprintf(fp_me_out, "%s ", name);
      fprintf(fp_mo_out, "%s ", name);
      if (metododireto==0 && print_input_pdfs==1)
	{
	  fprintf (fp_inppdfs_me, "%s ", name); fprintf (fp_inppdfs_mo, "%s ", name);
	}
      
      /* Checking null input values */
      if ( (teff <= -90.0 || steff <= -90.0) ||
	   (use_lg == 1.0 && (lg <= -90.0 || slg <= -90.0)) ||
	   (use_dnu == 1.0 && (dnu <= -90.0 || sdnu <= -90.0)) ||
	   (kind_tracks == 3 && use_pspacing == 1.0 && (pspac <= -90.0 || spspac <= -90.0)) || // only for seismic tracks: MESA
	   (use_numax == 1.0 && (numax <= -90.0 || snumax <= -90.0)) ||
	   (use_lum == 1.0 && (lum <= -90.0 || elum <= -90.0)) )
	{
	  printf("\nMissing defined INPUT parameters! :(\n");
	  for (i=0; i<48+(n_fil*5); i++)
	    fprintf(fp_mo_out, "%.4f ", ANNUL);
	  for (i=0; i<45+(n_fil*5); i++)	
	    fprintf(fp_me_out, "%.4f ", ANNUL);
	  fprintf (fp_mo_out, "0 null\n");
	  fprintf (fp_me_out, "0 null\n");
	  if (metododireto==0 && print_input_pdfs==1)
	    {
	      int nip=10; //Teff+MH
	      if (use_dnu>0.0) nip+=5;
	      if (kind_tracks==3 && use_pspacing>0.0) nip+=5;
	      if (use_numax>0.0) nip+=5;
	      if (use_lum>0.0) nip+=5;
	      for (i=0; i<nip;i++)
		{
		  fprintf (fp_inppdfs_me, "%.4f ", ANNUL); fprintf (fp_inppdfs_mo, "%.4f ", ANNUL);
		}
	      fprintf (fp_inppdfs_me, "\n"); fprintf (fp_inppdfs_mo, "\n");
	    }
	  continue;
	}
 
      /* Setting mu_min and max (HIPPARCOS-NADEGE = 1; APOKASC-CoRoGEE = 0) */
      bright_star = 0;
      for (ic=1; ic<n_fil+1; ic++)
	if (magobs.fil[ic] <= 6.0 && magobs.fil[ic] >= -90.0)
	  {
	    bright_star = 1;
	    break;
	  }
      if (bright_star == 1)
	{
	  mu0_min = -2.0;  mu0_max = 12.0;
	}
      else 
	{
	  mu0_min = 4.0;  mu0_max = 18.0; //TSR mu0_min=2.0 -> 4.0
	}
      
      /* Setting fitting extinction: use input Av = 0, fit Av = 1 */    
      fit_av=1;
      avx = 0.0; //avr
      savx = 0.0; //savr
      (fit_av == 1) ? printf ("\nFITTING EXTINCTION") : printf ("\nUSING INPUT EXTINCTION");
      
      /* Setting DELTA_TEFF */
      teff += delta_teff;
      if (delta_teff > 0.0)
	printf ("\nDelta_teff = %.2f, teff = %.2f K -> %.2f K",
		delta_teff, teff, teff+delta_teff);
    
      /* Printing null magnitudes (only to have control)  */
      printf("\nNULL magnitudes: ");
      for (ic=1; ic<n_fil+1; ic++)
	if (magobs.fil[ic] <= -90.0)
	  printf("%s ", nomi_mag[ic]);
      
      /* probability distributions to be looked for: */
      //set_pdf(&pdf, flag_log, xmin, xmax, dx);
      set_pdf(&pdf_age, 1, log10(agemin), log10(agemax), 0.01); //TSR: MODIFIQUEI DELTA_AGE 0.02->0.01
      set_pdf(&pdf_m, 1, log10(0.5), log10(6.0), 0.002); //TSR: MODIFIQUEI DELTA_MASS 0.005->0.002 29/07/2016
      set_pdf(&pdf_g, 1, -2.0, 5.0, 0.005); 
      set_pdf(&pdf_r, 1, -1.0, 2.0, 0.004); 
      set_pdf(&pdf_rho, 1, log10(1e-5), log10(2.0), 0.0025);
      for (ic=0; ic<n_fil+1; ic++) 
	{
	  set_pdf(&pdf_mag[ic], 0, -12.0, 10.0, 0.020); // intervalo de mag grande!
	  set_pdf(&pdf_mu[ic], 0, mu0_min, mu0_max, 0.005);
	}
      set_pdf(&pdf_av, 0, -1.5, 7.5, 0.01); //TSR: modifiquei limite superior de Av 1.0->3.0->7.5 e inferior -0.5->-1.5
      set_pdf(&pdf_mu0, 0, mu0_min, mu0_max, 0.005);
      //Input_pdfs: need to be improved!
      if (metododireto==0 && print_input_pdfs==1) 
	{
	  set_pdf(&pdf_logte, 1, 3.2, 4.0, 0.001); //TSR:MODIFIQUEI DELTA_TEFF 0.002->0.001 04/05/2015
	  set_pdf(&pdf_fehv, 0, -3.0, 1.0, 0.01);
	  if (use_dnu>0.0) set_pdf(&pdf_deltanu, 0, (dnu-10.0*sdnu <0) ? 0.0 : (int) (dnu-10.0*sdnu),
				   (int) (dnu+10.0*sdnu), (sdnu <0.0) ? (int) ((sdnu*100)/100) : 0.2); //TSR: MODIFIQUEI DETLA_DNU 0.005 ->
	  if (kind_tracks==3 && use_pspacing>0.0)
	    set_pdf(&pdf_pspacing, 0, (pspac-10.0*spspac <0) ? 0.0 : (int) (pspac-10.0*spspac),
		    (int) (pspac+10.0*pspac), (int) (0.5*spspac)); 
	  if (use_numax>0.0) set_pdf(&pdf_numax0, 0, (numax-10.0*snumax <0) ? 0.0 : (int) (numax-10.0*snumax),
				     (int) (numax+10.0*snumax), (int) (0.2*snumax));
	  // printf("\nDnu min %f, max %f delta %f \n", pdf_deltanu.xmin, pdf_deltanu.xmax, pdf_deltanu.dx);
	  // printf("\nnumax min %f, max %f delta %f \n", pdf_numax0.xmin, pdf_numax0.xmax, pdf_numax0.dx);
	  // printf("\npspacing min %f, max %f delta %f \n", pdf_pspacing.xmin, pdf_pspacing.xmax, pdf_pspacing.dx);
	  if (use_lum>0.0) set_pdf(&pdf_loglum, 1, -1.0, 4.0, 0.001);
	}
      //JPDFS
      if (print_jpdfs==1)
	{
	  set_pdf2(&jpdf_age_m, &pdf_age, &pdf_m);
	  set_pdf2(&jpdf_age_g, &pdf_age, &pdf_g);
	  set_pdf2(&jpdf_age_r, &pdf_age, &pdf_r);
	  set_pdf2(&jpdf_m_g, &pdf_m, &pdf_g);
	  set_pdf2(&jpdf_m_r, &pdf_m, &pdf_r);
	  set_pdf2(&jpdf_r_g, &pdf_r, &pdf_g);
	}
      //Extinction table
      if (ext_coef_table==1)
	{
	  alamb_av_teff (teff, &alamb_av[0]) ;
	  if (print_pdfs==1)
	    {
	      strcpycat (temp, PDFSDIR, "AlambAv/");
	      strcat (temp, "alamb_av_");
	      strcat (temp, name);
	      fp_avlamb = fopen (temp, "w") ;
	      fprintf (fp_avlamb, "#Id " ) ;
	      for (ic=1; ic<n_fil+1; ic++)
		fprintf (fp_avlamb, "%s ", nomi_mag[ic]) ;
	      
	      fprintf (fp_avlamb, "\n%s ", name);
	      for (ic=1; ic<n_fil+1; ic++)
		fprintf (fp_avlamb, "%f ", alamb_av[ic]) ;
	      fclose(fp_avlamb) ;
	    }
	}
      
#ifdef SCRIVE_ISOC
      FILE *fp_isoc, *fp_logl, *fp_jpfloglte;
      PDF pdf_logl;
      PDF2 jpdf_logl_logte;
      float isocage=-99.9, isoczeta=-99.9;
      strcpy(temp, "isocs_int_");
      strcat (temp, name);
      fp_isoc = fopen(temp, "a");
      //fp_logl = fopen("pdflogl.test", "w");      
      set_pdf(&pdf_logl, 1, -1.0, 4.0, 0.002);
      set_pdf2(&jpdf_logl_logte, &pdf_logl, &pdf_logte);
#endif

      if (teff<=0.0) continue;
      
      lteff = log10(teff) ;
      slteff = steff / (teff * log(10.0)) ;
      
#ifdef ANDREAOUTPUT
      // Andrea 
      // opens file where PDF will be written (1 file per star)
      char outputfile2[200];
      strcpy (outputfile2, name );
      strcat(outputfile2, ".pdf.txt");
      fp_out_tmp=fopen(outputfile2,"w");
      printf("doing %s\n",outputfile2);
      probtot = 0.0;
      for (im=0; im<pdf_m.n; im++) 
	for (ir=0; ir<pdf_r.n; ir++)
	  for (j=0; j<pdf_age.n; j++)
	    prob[im][ir][j] = 0.0 ;
#endif      
      
      /* loop over feh */
      if (mh>-90.0) // if there's an observed metallicity, use it 
	{
	  feh=mh;
	  sfeh=smh;
	  delta_feh = 0.01 ; /* 0.01 dex default resolution */
	  feh_min = feh-((int)(NSIGMA*sfeh/delta_feh))*delta_feh ; 
	  feh_max = feh+((int)(NSIGMA*sfeh/delta_feh))*delta_feh ;
	}
      else // use a very broad distribution as prior
 	{
	  feh=-0.5;
	  sfeh=1.0;
	  delta_feh = 0.01 ; /* 0.01 dex default resolution */
	  feh_min = -2.3 ; 
	  feh_max = +0.3 ;
	}

      //printf("\nZ %f", zfunc(mh));
      //printf("\nmh %f", mhfunc(0.016));
      if (metododireto==0) // metodo Bayesiano
	{
	  //////// main loop over [Fe/H] and log(age): //////////////
	  for (fehv=feh_min; fehv<=feh_max; fehv += delta_feh )  
	    {
	      // printf("%f\n",fehv);
	      // simply assumes constant weight for different [Fe/H],
	      // OR use measured values (and sigma) from APOGEE
	      fact_feh = (sfeh==0.0) ?  1.0 : gaussiana (1.0, sfeh, feh-fehv) ;
	      zeta = zfunc(fehv);//0.019 * pow(10.0,fehv) ;
	      /* acha probabilidade baseada em menor distancia */
	      for (lage=pdf_age.xmin, j=0; j<pdf_age.n; lage+=pdf_age.dx, j++) 
		{
		  age = pow (10.0, lage) ;
		  //printf(" %f %g\n",fehv, age);
		  if (age>agemax) continue ;
		  // fact_age is simply Delta(age) in years, for each age bin
		  fact_age = pow (10.0, lage+pdf_age.dx) - age ;
		  fa_isocrona_completa (zeta, age, isoc);
		  mass_ini = isoc->mini[0] ;
		  mass_act = isoc->star[0].mass ;
		  logl = isoc->star[0].logl ;
		  logte = isoc->star[0].logte ;
		  if (kind_tracks==3) // only for seismic tracks: MESA
		    {
		      //deltanu = isoc->star[0].deltanu ;
		      deltanu = isoc->star[0].ddeltanu * fdeltanu(dnu_sun_obs, mass_act, pow(10.0, isoc->star[0].logr)); // interpolating Dnu/DnuSR
		      pspacing = isoc->star[0].pspacing ;
		    }
		  
		  for (i=0; i<isoc->npun-1; i++)
		    {
		      // printf ("%d\n",i);
		      /* n_int = HR_RES * (1 + (int)(fabs(isoc->logl[i+1]- */
		      /* 			     isoc->logl[i])/d_logl) ); */
		      n_int = HR_RES ;
		      delta_mass_ini = (isoc->mini[i+1] - isoc->mini[i]) / 
			(float)n_int ;
		      delta_mass_act = (isoc->star[i+1].mass - isoc->star[i].mass) / 
			(float)n_int ;
		      delta_logl = (isoc->star[i+1].logl - isoc->star[i].logl) / 
			(float)n_int;
		      delta_logte = (isoc->star[i+1].logte - isoc->star[i].logte) / 
			(float)n_int;
		      if (kind_tracks==3) // only for seismic tracks: MESA
			{
			  // delta_dnu = (isoc->star[i+1].deltanu - isoc->star[i].deltanu) / 
			  // (float)n_int; 
			  delta_dnu = (isoc->star[i+1].ddeltanu * fdeltanu(dnu_sun_obs, isoc->star[i+1].mass, pow(10.0, isoc->star[i+1].logr)) - 
				       isoc->star[i].ddeltanu * fdeltanu(dnu_sun_obs, isoc->star[i].mass, pow(10.0, isoc->star[i].logr)) ) / (float)n_int; 
			  delta_pspacing = (isoc->star[i+1].pspacing - isoc->star[i].pspacing) / 
			    (float)n_int;
			}
		      
		      for (ii=0; ii<n_int; ii++)
			{
      			  logg = loggrav (mass_act, logl, logte) ;
			  // trasforma_in_magnitu (zeta, logl, logg, logte, 
			  //			0.0, -1, 0.0,
			  //			CO_DEFAULT, 0.0, &magni);
			  star_magnitu (zeta, logl, logg, logte,  //TSR ADDED 26.01.2016
					&isoc->star[i], //WARNING: CHECK THIS!!!!!!!
					0.0, &magni);
			  mv0 = magni.fil[3] ;
			  teff0 = pow(10.0, logte) ;
			  logr = flogradius(logl, logte) ; 
			  r0 = pow(10.0, logr) ;
			  lrho = flogrho(mass_act, logr);
			  if (kind_tracks!=3) // only for NON-seismic tracks (1 and 2)
			    deltanu = fdeltanu(dnu_sun_obs, mass_act, r0) ;		  
			  numax0 = fnumax(mass_act, r0, teff0) ;
			  
			  if ((use_teff*fabs(lteff-logte)<NSIGMA*slteff) && 
			      (use_lg*fabs(lg-logg)<fabs(NSIGMA*slg)) &&
			      (use_dnu*fabs(dnu-deltanu)<fabs(NSIGMA*sdnu)) &&
			      // only for seismic tracks: MESA
			      (use_pspacing*fabs(pspac-pspacing)<fabs(NSIGMA*spspac)) &&
			      (use_numax*fabs(numax-numax0)<fabs(NSIGMA*snumax)) &&
			      // weigths only stars close to obs. Teff,numax,dnu
			      // this is the place to add new observables, e.g. period spacing:
			      ( (evstate==1 && isoc->star[i].label<=LABEL_RGB) || //on RGB
				(evstate==2 && isoc->star[i].label>=LABEL_CHEB && isoc->star[i].label<LABEL_EAGB) || //on CHeB
				(evstate==3 && isoc->star[i].label<=LABEL_RGB && isoc->star[i].label>=LABEL_EAGB) || //on RGB or AGB -- TSR added on 08/08/2016
				(evstate!=1 && evstate!=2 && evstate!=3) //non identified phase
				)
			      )
			    {
			      // this is the weight for quantities that vary along the isochrones:
			      //printf ("\ni+1 %.4f, i %.4f", isoc->mass_act[i+1], isoc->mass_act[i]);
			      //gaussiana (1.0, observational error, observational-model)
			      weight = 
				((use_teff>0.0) ? use_teff * gaussiana (1.0, slteff, lteff-logte) : 1.0 ) *
				((use_lg>0.0) ? use_lg * gaussiana (1.0, slg, lg-logg) : 1.0 ) *				
				((use_lum>0.0) ? use_lum * gaussiana (1.0, slog_lum, log_lum-logl) : 1.0 ) *
				((use_dnu>0.0) ? use_dnu * gaussiana (1.0, sdnu, dnu-deltanu) : 1.0 ) *
				((use_lum>0.0) ? use_lum * gaussiana (1.0, 0.9, 44.2-r0) : 1.0 ) *
				//printf ("using inteferometric constraints 44.2R +- 0.9")

				//Added line to include inteferometric radius constraint if wanted: Thomas North 19/9/17
				// only for seismic tracks: MESA
				((kind_tracks==3 && use_pspacing>0.0) ? use_pspacing * gaussiana (1.0, spspac, pspac-pspacing) : 1.0 ) *
				((use_numax>0.0) ? use_numax * gaussiana (1.0, snumax, numax-numax0) : 1.0 ) ;
			      // printf ("weight %f\n", weight) ;
			      // this is the weight for quantities that do NOT vary along isochrones:
			      nstelle = (fact_age/1.0e9) * fact_feh *
				imf_n (mass_ini, mass_ini+delta_mass_ini) ;
			      if (nstelle>0.0) {
				add2pdf (&pdf_age, lage, nstelle*weight);
				add2pdf (&pdf_m, log10(mass_act), nstelle*weight);
				add2pdf (&pdf_g, logg, nstelle*weight);
				add2pdf (&pdf_r, logr, nstelle*weight);
				add2pdf (&pdf_rho, lrho, nstelle*weight);
				for (ic=0; ic<n_fil+1; ic++)
				  add2pdf (&pdf_mag[ic], magni.fil[ic], nstelle*weight);
				if (print_input_pdfs==1)
				  {
				    add2pdf (&pdf_logte, logte, nstelle*weight);
				    add2pdf (&pdf_fehv, fehv, nstelle*weight);
				    if (use_dnu>0.0) add2pdf (&pdf_deltanu, deltanu, nstelle*weight);
				    if (kind_tracks==3 && use_pspacing>0.0)
				      add2pdf (&pdf_pspacing, pspacing, nstelle*weight);
				    if (use_numax>0.0) add2pdf (&pdf_numax0, numax0, nstelle*weight);
				    if (use_lum>0.0) add2pdf (&pdf_loglum, logl, nstelle*weight);
				  }
#ifdef SCRIVE_ISOC
				add2pdf (&pdf_logl, logl, nstelle*weight);
				add2pdf2 (&jpdf_logl_logte, logl, logte, nstelle*weight);
				if (isocage!=age || isoczeta!=zeta) scrive_isocrona (fp_isoc, isoc, 0, 10.0, 1.0);
				isocage=age;
				isoczeta=zeta;
#endif
				
				if (print_jpdfs==1)
				  {
				    add2pdf2 (&jpdf_age_m, lage, log10(mass_act), nstelle*weight);
				    add2pdf2 (&jpdf_age_g, lage, logg, nstelle*weight);
				    add2pdf2 (&jpdf_age_r, lage, logr, nstelle*weight);
				    add2pdf2 (&jpdf_m_g, log10(mass_act), logg, nstelle*weight);
				    add2pdf2 (&jpdf_m_r, log10(mass_act), logr, nstelle*weight);
				    add2pdf2 (&jpdf_r_g, logr, logg, nstelle*weight);
				  }
				
#ifdef ANDREAOUTPUT		
				im = (int)( (log10(mass_act)-lm_ini)/delta_lm ) ;
				ir = (int)( (logr-lr_ini)/delta_lr ) ;
				if (im>=0 && im<pdf_m.n &&
				    ir>=0 && ir<pdf_r.n &&
				    j>=0 && j<pdf_age.n) 
				  {
				    prob[im][ir][j] += nstelle * weight ;
				    probtot +=  nstelle * weight ;
				  }
#endif
			      }
			    }
			  
			  mass_ini += delta_mass_ini ;
			  mass_act += delta_mass_act ;
			  logl += delta_logl ;
			  logte += delta_logte ;
			  if (kind_tracks==3) // only for seismic tracks: MESA
			    {
			      deltanu += delta_dnu ;
			      pspacing += delta_pspacing ;
			    }		  
			}
		    }
		}
	    }
	}
      else // metodo direto: sets all pdfs as gaussians
	{
	  STAR stnull ; //TSR ADDED 26.01.2016
	  starnull(stnull);  //TSR ADDED 26.01.2016
	  zeta = zfunc(feh);//0.019 * pow(10.0,feh) ;
	  metdir (dnu, sdnu, numax, snumax, teff, steff,
		  &logrhodir, &elogrhodir, &massdir, &emassdir, 
		  &raiodir, &eraiodir, &loggdir, &eloggdir);
	  /* printf ("\n m = %.3f +- %.3f Ms\n r = %.3f +- %.3f Rs\nlog g =  %.3f +- %.3f\nlog rho = %.3f +- %.3f\n",
	     massdir, emassdir, raiodir, eraiodir,
	     loggdir, eloggdir, logrhodir, elogrhodir);
	  */
	  logmassdir = log10(massdir);
	  elogmassdir = emassdir/(massdir*log(10.0)) ;
	  lograiodir = log10(raiodir);
	  elograiodir = eraiodir/(raiodir*log(10.0)) ;
	  set_pdf2gaussian (&pdf_age, log10(0.5*(agemin+agemax)), 1.0); // meaningless
	  set_pdf2gaussian (&pdf_m, (float)logmassdir, (float)elogmassdir); 
	  set_pdf2gaussian (&pdf_g, (float)loggdir, (float)eloggdir); 
	  set_pdf2gaussian (&pdf_r, (float)lograiodir, (float)elograiodir); 
	  set_pdf2gaussian (&pdf_rho, (float)logrhodir, (float)elogrhodir);
	  logl = 2.0*(float)lograiodir // log10(R/Rsun)^2 
	    + 4.0*((float)lteff-logte_sun) ; // log10(Teff/Teff_sun)^4
	  /* trasforma_in_magnitu (zeta, logl, (float)loggdir, (float)lteff, 
	     0.0, -1, 0.0,
	     CO_DEFAULT, 0.0, &magni);
	     trasforma_in_magnitu (zeta, logl, (float)loggdir, 
	     (float)(lteff+slteff), 
	     0.0, -1, 0.0,
	     CO_DEFAULT, 0.0, &deltamag); */
	  star_magnitu (zeta, logl, (float)loggdir, (float)lteff,  //TSR ADDED 26.01.2016
			&stnull,
			0.0, &magni);
	  star_magnitu (zeta, logl, (float)loggdir,  //TSR ADDED 26.01.2016
			(float)(lteff+slteff), 
			&stnull,
			0.0, &deltamag);
	  for (ic=0; ic<n_fil+1; ic++)
	    {
	      emag = sqrt(
			  // error in magnitude due to errors in radius and Teff only:
			  pow(5.0*(float)elograiodir, 2.0) +
			  pow(10.0*(float)slteff, 2.0) +
			  // add error due to BC(Teff)
			  pow(magni.fil[ic]-deltamag.fil[ic], 2.0)
			  );
	      set_pdf2gaussian (&pdf_mag[ic], magni.fil[ic], emag); // PDF for absolute mag
	    }
	}
#ifdef SCRIVE_ISOC
      fclose(fp_isoc);
#endif 
      printf ("\n\n%s", name);
    
      /* PDF logAGE: median, mode and CI of 68%, 95% */
      //P(logt)dlogt = P(t)dt => P(t)=P(logt)(dlogt/dt)
      //t=10^logt => P(t) = P(logt) / (ln10 * t) 
      if (metododireto==0)
	{
	  computeCI_pdf (&pdf_age);
	  print_PDFsummary (&pdf_age, "age", "Gyr");
	  print_outputfiles (fp_mo_out, fp_me_out, &pdf_age, 0);
	  
	}
      else
	{
	  fprintf(fp_me_out, "%.4f %.4f %.4f %.4f %.4f ", ANNUL, ANNUL,
		  ANNUL, ANNUL, ANNUL);
	  fprintf(fp_mo_out, "%.4f %.4f %.4f %.4f %.4f ", ANNUL, ANNUL,
		  ANNUL, ANNUL, ANNUL);
	}
      
      /* PDF logMASS: median, mode and CI of 68%, 95% */
      if (metododireto==0)
	{
	  computeCI_pdf (&pdf_m);
	  print_PDFsummary (&pdf_m, "m", "Ms");
	  print_outputfiles (fp_mo_out, fp_me_out, &pdf_m, 0);
	}
      else
	{
	  printf("\n m  = %.2f Ms, CI(68%%): %.2f - %.2f, CI(95%%): %.2f - %.2f", 
		 massdir, (massdir-emassdir), (massdir+emassdir),
		 (massdir-(2*emassdir)), (massdir+(2*emassdir)));
	  fprintf(fp_me_out, "%.4f %.4f %.4f %.4f %.4f ", massdir,
		  (massdir-emassdir), (massdir+emassdir),
		  (massdir-(2*emassdir)), (massdir+(2*emassdir)));
	  fprintf(fp_mo_out, "%.4f %.4f %.4f %.4f %.4f ", massdir,
		  (massdir-emassdir), (massdir+emassdir),
		  (massdir-(2*emassdir)), (massdir+(2*emassdir)));
	}
        
      /* PDF logg: median, mode and CI of 68%, 95% */
      if (metododireto==0)
	{
	  computeCI_pdf (&pdf_g);
	  print_PDFsummary (&pdf_g, "g", "0");
	  print_outputfiles (fp_mo_out, fp_me_out, &pdf_g, 1);
	}
      else 
	{
	  printf("\n logg = %.2f, CI(68%%): %.2f - %.2f, CI(95%%): %.2f - %.2f",
		 loggdir, (loggdir-eloggdir), (loggdir+eloggdir),
		 (loggdir-(2.0*eloggdir)), (loggdir+(2.0*eloggdir)));
	  fprintf(fp_me_out, "%.4f %.4f %.4f %.4f %.4f ", loggdir,
		  (loggdir-eloggdir), (loggdir+eloggdir),
		  (loggdir-(2.0*eloggdir)), (loggdir+(2.0*eloggdir)));
	  fprintf(fp_mo_out, "%.4f %.4f %.4f %.4f %.4f ",loggdir,
		  (loggdir-eloggdir), (loggdir+eloggdir),
		  (loggdir-(2.0*eloggdir)), (loggdir+(2.0*eloggdir)));
	}
      
      /* PDF logRADIUS: median, mode and CI of 68%, 95% */
      if (metododireto==0)
	{
	  computeCI_pdf (&pdf_r);
	  print_PDFsummary (&pdf_r, "r", "Rs");
	  print_outputfiles (fp_mo_out, fp_me_out, &pdf_r, 0);
	}
      else 
	{	 
	  printf("\n r = %.2f Rs, CI(68%%): %.2f - %.2f, CI(95%%): %.2f - %.2f", 
		 raiodir, (raiodir-eraiodir), (raiodir+eraiodir),
		 (raiodir-(2.0*eraiodir)), (raiodir+(2.0*eraiodir)));
	  fprintf(fp_me_out, "%.4f %.4f %.4f %.4f %.4f ", raiodir,
		  (raiodir-eraiodir), (raiodir+eraiodir),
		  (raiodir-(2.0*eraiodir)), (raiodir+(2.0*eraiodir)));
	  fprintf(fp_mo_out, "%.4f %.4f %.4f %.4f %.4f " ,raiodir,
		  (raiodir-eraiodir), (raiodir+eraiodir),
		  (raiodir-(2.0*eraiodir)), (raiodir+(2.0*eraiodir)));
	}
             
      /* PDF logRHO: median, mode and CI of 68%, 95% */
      if (metododireto == 0)
	{
	  computeCI_pdf (&pdf_rho);
	  print_PDFsummary (&pdf_rho, "rho", "0");
	  print_outputfiles (fp_mo_out, fp_me_out, &pdf_rho, 1);
	}
      else
	{
	  printf("\n logrho = %.2f, CI(68%%): %.2f - %.2f, CI(95%%): %.2f - %.2f",
		 logrhodir, (logrhodir-elogrhodir), (logrhodir+elogrhodir),
		 (logrhodir-(2.0*elogrhodir)), (logrhodir+(2.0*elogrhodir)));
	  fprintf(fp_me_out, "%.4f %.4f %.4f %.4f %.4f ", logrhodir,
		  (logrhodir-elogrhodir), (logrhodir+elogrhodir),
		  (logrhodir-(2.0*elogrhodir)), (logrhodir+(2.0*elogrhodir)));
	  fprintf(fp_mo_out, "%.4f %.4f %.4f %.4f %.4f ", logrhodir,
		  (logrhodir-elogrhodir), (logrhodir+elogrhodir),
		  (logrhodir-(2.0*elogrhodir)), (logrhodir+(2.0*elogrhodir)));
	}

      /* Printing PDFs */
      if (metododireto==0 && print_pdfs==1)
	{
	  print_PDFfile(fp_age, "PDFage/", "pdf_age_", name, "#logage PDFlogage age(Gyr) PDF\n",
			&pdf_age, 1.0e-8, 0);
	  print_PDFfile(fp_mass, "PDFm/", "pdf_m_", name, "#logm PDFlogm m(Ms) PDF\n",
			&pdf_m, 1.0e-9, 0);
	  print_PDFfile(fp_logg, "PDFlogg/", "pdf_logg_", name, "#logg PDF\n",
			&pdf_g, 1.0e-9, 1);
	  print_PDFfile (fp_r, "PDFr/", "pdf_r_", name, "#logr PDFlogr r(Rs) PDF\n",
			 &pdf_r, 1.0e-8, 0);
	  print_PDFfile (fp_rho, "PDFlrho/", "pdf_lrho_", name, "#logrho PDF\n",
			 &pdf_rho, 1.0e-7, 1);
#ifdef SCRIVE_ISOC
	  normalizes_pdf(&pdf_logl);
	  print_PDFfile (fp_logl, "PDFlogL/", "pdf_logL_", name, "#logL_Lsun PDF\n",
			 &pdf_logl, 1.0e-7, 1);
#endif
	}
      /* Printing INPUT PDFs */
      if (metododireto==0 && print_input_pdfs==1)
	{
	  computeCI_pdf (&pdf_logte);
	  print_outputfiles (fp_inppdfs_mo, fp_inppdfs_me, &pdf_logte, 0);
	  print_PDFfile (fp_logte, "PDFlogte/", "pdf_logte_", name, "#logTeff PDF\n",
			 &pdf_logte, 1.0e-9, 1);
	  
	  computeCI_pdf (&pdf_fehv);
	  print_outputfiles (fp_inppdfs_mo, fp_inppdfs_me, &pdf_fehv, 0);
	  print_PDFfile (fp_fehv, "PDFmh/", "pdf_mh_", name, "#[M/H] PDF\n",
			 &pdf_fehv, 1.0e-9, 1);
	  if (use_dnu>0.0)
	    {
	      computeCI_pdf (&pdf_deltanu);
	      print_outputfiles (fp_inppdfs_mo, fp_inppdfs_me, &pdf_deltanu, 0);
	      print_PDFfile (fp_deltanu, "PDFDnu/", "pdf_dnu_", name, "#Dnu PDF\n",
			     &pdf_deltanu, 1.0e-9, 1);
	    }
	  if (kind_tracks==3 && use_pspacing>0.0)
	    {
	      computeCI_pdf (&pdf_pspacing);
	      print_outputfiles (fp_inppdfs_mo, fp_inppdfs_me, &pdf_pspacing, 0);
	      print_PDFfile (fp_pspacing, "PDFDPi1/", "pdf_dpi1_", name, "#DPi1 PDF\n",
			     &pdf_pspacing, 1.0e-9, 1);
	    }
	  if (use_numax>0.0)
	    {
	      computeCI_pdf (&pdf_numax0);
	      print_outputfiles (fp_inppdfs_mo, fp_inppdfs_me, &pdf_numax0, 0);
	      print_PDFfile (fp_numax0, "PDFnumax/", "pdf_numax_", name, "#numax PDF\n",
			     &pdf_numax0, 1.0e-9, 1);
	    }
          if (use_lum>0.0)
            {
	      computeCI_pdf (&pdf_loglum);
	      print_outputfiles (fp_inppdfs_mo, fp_inppdfs_me, &pdf_loglum, 1);
              print_PDFfile (fp_loglum, "PDFlogLum/", "pdf_logLum_", name, "#logL_Lsun PDF\n",
                             &pdf_loglum, 1.0e-7, 1);
            }
	  fprintf (fp_inppdfs_me, "\n"); fprintf (fp_inppdfs_mo, "\n");
	}      
      /* Printing JPDFs */
      if (metododireto==0 && print_jpdfs==1)
	{
	  normalizes_pdf2(&jpdf_age_m);
	  print_PDF2file (fp_jpdflogagem, "JPDF_age_m/", "jpdf_age_m_", name, "#logage logm JPDF\n",
			  &jpdf_age_m, 1.0e-10);
	  normalizes_pdf2(&jpdf_age_g);
	  print_PDF2file (fp_jpdflogageg, "JPDF_age_g/", "jpdf_age_g_", name, "#logage logg JPDF\n",
			  &jpdf_age_g, 1.0e-10);
	  normalizes_pdf2(&jpdf_age_r);
	  print_PDF2file (fp_jpdflogager, "JPDF_age_r/", "jpdf_age_r_", name, "#logage logr JPDF\n",
			  &jpdf_age_r, 1.0e-10);
	  normalizes_pdf2(&jpdf_m_g);
	  print_PDF2file (fp_jpdflogmg, "JPDF_m_g/", "jpdf_m_g_", name, "#logm logg JPDF\n",
			  &jpdf_m_g, 1.0e-10);
       	  normalizes_pdf2(&jpdf_m_r);
	  print_PDF2file (fp_jpdflogmr, "JPDF_m_r/", "jpdf_m_r_", name, "#logm logr JPDF\n",
			  &jpdf_m_r, 1.0e-10);
      	  normalizes_pdf2(&jpdf_r_g);
	  print_PDF2file (fp_jpdflogrg, "JPDF_r_g/", "jpdf_r_g_", name, "#logr logg JPDF\n",
			  &jpdf_r_g, 1.0e-10);
#ifdef SCRIVE_ISOC
	  normalizes_pdf2(&jpdf_logl_logte);
	  print_PDF2file (fp_jpfloglte, "JPDF_logL_logte/", "jpdf_logL_logte_", name, "#logL_Lsun logte JPDF\n",
			 &jpdf_logl_logte, 1.0e-10);
#endif
	}
      
#ifdef thaise_teste
      float dvol, volum, dvoll, voluml;
      for (im=1; im<pdf_m.n-1; im++) 
	for (ir=1; ir<pdf_r.n-1; ir++)
	  {
	    dvol = 1.0;
	    dvol *= (pow(10.0, pdf_m.xmin+(float)(im+1)*pdf_m.dx)-
		     pow(10.0, pdf_m.xmin+(float)(im)*pdf_m.dx))/2.0+
	      (pow(10.0, pdf_m.xmin+(float)(im)*pdf_m.dx)-
	       pow(10.0, pdf_m.xmin+(float)(im-1)*pdf_m.dx))/2.0;
	    
	    dvol *= (pow(10.0, pdf_r.xmin+(float)(ir+1)*pdf_r.dx)-
		     pow(10.0, pdf_r.xmin+(float)(ir)*pdf_r.dx))/2.0+
	      (pow(10.0, pdf_r.xmin+(float)(ir)*pdf_r.dx)-
	       pow(10.0, pdf_r.xmin+(float)(ir-1)*pdf_r.dx))/2.0;
	    dvol *=  (jpdf_logm_logr[im][ir]/probtotmr)*
	      (log(10.0)*log(10.0)*pow(10.0, pdf_m.xmin+(float)(im)*pdf_m.dx)*
	       pow(10.0, pdf_r.xmin+(float)(ir)*pdf_r.dx));
	    volum += dvol;
	    
	    dvoll = 1.0;
	    dvoll *= abs((pdf_m.xmin+(float)(im+1)*pdf_m.dx)-  printf("\nfile: %s", temp);

			(pdf_m.xmin+(float)(im)*pdf_m.dx))/2.0+
	      abs((pdf_m.xmin+(float)(im)*pdf_m.dx)-
		  (pdf_m.xmin+(float)(im-1)*pdf_m.dx))/2.0;
	    dvoll *= abs((pdf_r.xmin+(float)(ir+1)*pdf_r.dx)-
			 (pdf_r.xmin+(float)(ir)*pdf_r.dx))/2.0+
	      abs((pdf_r.xmin+(float)(ir)*pdf_r.dx)-
	       (pdf_r.xmin+(float)(ir-1)*pdf_r.dx))/2.0;
	    dvoll *= (jpdf_logm_logr[im][ir]/probtotmr);
	    voluml += dvoll;
	  }

      printf("\nvolume: %f\n", volum);
      printf("\nvolumel: %f\n", voluml);

      for (im=0; im<pdf_m.n; im++) 
	for (ir=0; ir<pdf_r.n; ir++)
	  {
	    // if (prob[im][ir][j]/probtot > 1.0e-10)
	    fprintf(fp_teste, "%d %d %g %g %g\n", im, ir,  
		    pdf_m.xmin+(float)(im)*pdf_m.dx,
		    pdf_r.xmin+(float)(ir)*pdf_r.dx,
		    jpdf_logm_logr[im][ir]/probtotmr,
		    pow(10.0, pdf_m.xmin+(float)(im)*pdf_m.dx),
		    pow(10.0, pdf_r.xmin+(float)(ir)*pdf_r.dx),
		    (jpdf_logm_logr[im][ir]/probtotmr)*
		    (log(10.0)*log(10.0)*pow(10.0, pdf_m.xmin+(float)(im)*pdf_m.dx)*
		     pow(10.0, pdf_r.xmin+(float)(ir)*pdf_r.dx)));
	  }

      fclose(fp_teste);
#endif

#ifdef ANDREAOUTPUT
      // prints big matrix with M,R,age probability:
      char outputfile3[200];
      strcpy (outputfile3, name );
      strcat(outputfile3, ".pppdf.txt");
      fppp_out=fopen(outputfile3,"w");
      printf("doing %s\n",outputfile3);
      fprintf(fppp_out, "# Total 3D Prob. for %s = %.5g\n", name, probtot); 
      fprintf(fppp_out, "# logM min, max, delta, n: %.4f %.4f %.4f %d\n", 
	      pdf_m.xmin, pdf_m.xmin+(float)(pdf_m.n)*pdf_m.dx, pdf_m.dx, pdf_m.n);       
      fprintf(fppp_out, "# logR min, max, delta, n: %.4f %.4f %.4f %d\n", 
	      pdf_r.xmin, pdf_r.xmin+(float)(pdf_r.n)*pdf_r.dx, pdf_r.dx, pdf_r.n);       
      fprintf(fppp_out, "# logt min, max, delta, n: %.4f %.4f %.4f %d\n", 
	      pdf_age.xmin, pdf_age.xmin+(float)(pdf_age.n)*pdf_age.dx, pdf_age.dx, pdf_age.n); 
      for (im=0; im<pdf_m.n; im++) 
	for (ir=0; ir<pdf_r.n; ir++)
	  for (j=0; j<pdf_age.n; j++)
	    {
	      if (prob[im][ir][j]/probtot > 1.0e-10)
		fprintf(fppp_out, "%d %d %d %.3g %.3g %.3g %.5g\n", im, ir, j,  
			pow(10.0,pdf_m.xmin+(float)(im)*pdf_m.dx),
			pow(10.0,pdf_r.xmin+(float)(ir)*pdf_r.dx),
			pow(10.0,pdf_age.xmin+(float)(j)*pdf_age.dx),
			prob[im][ir][j]/probtot) ;
	    }
      fclose(fppp_out);
#endif

      strcpy(nomimag,"");
      if (print_jpdfs==1)
	{
	  strcpycat (temp, PDFSDIR, "JPDF_mag/");
	  strcat (temp, "jpdf_mag_");
	  strcat (temp, name);
	  fp_jpdfmag = fopen (temp, "w") ;
	  fprintf (fp_jpdfmag, "#fil magobs mag mu JPDF\n");
	}
      for (nic=0, ic=0; ic<n_fil+1; ic++)
	{	  
	  /* PDF MAG: median, mode and CI of 68%, 95% */
    	  computeCI_pdf (&pdf_mag[ic]);
	  if (print_pdfs==1) print_PDFmagfile (fp_mag, "PDFmag/", "pdf_mag_", name, "#fil mag PDF\n",
					       &pdf_mag[ic], 1.0e-6, ic);
	  print_PDFsummary (&pdf_mag[ic], nomi_mag[ic], "0");
	  print_outputfiles (fp_mo_out, fp_me_out, &pdf_mag[ic], 0); 
	  
	  if (ic>0 && magobs.fil[ic] > -90.0 && smagobs.fil[ic] > -90.0)
	    { 
	      /* Guassian distribution magobs */
	      set_pdf(&pdf_magobs[ic], 0, magobs.fil[ic]-5.0*smagobs.fil[ic],
		      magobs.fil[ic]+5.0*smagobs.fil[ic], 0.005);
	      for (i=0; i<pdf_magobs[ic].n; i++)
		pdf_magobs[ic].y[i] = gaussiana (1.0, smagobs.fil[ic],
						 -5.0*smagobs.fil[ic]+(float)(i)*pdf_magobs[ic].dx);
	      normalizes_pdf(&pdf_magobs[ic]);
	      if (print_pdfs==1) print_PDFmagfile (fp_magobs, "PDFmagobs/", "pdf_magobs_", name, "#fil magobs PDF\n",
						   &pdf_magobs[ic], 1.0e-6, ic);
	      
	      /* JPDF (mag, magobs)  */
	      for (ijpdf=0, i=0; i<pdf_mag[ic].n; i++) //mag
		for (j=0; j<pdf_magobs[ic].n; j++) //magobs 
		  {
		    if ((pdf_magobs[ic].y[j]*pdf_mag[ic].y[i])>=1.0e-6 && print_jpdfs==1)
		      fprintf(fp_jpdfmag, "%d %lf %lf %lf %lf\n", ic, 
			      (pdf_magobs[ic].xmin+(float)(j)*pdf_magobs[ic].dx), //magobs
			      (pdf_mag[ic].xmin+(float)(i)*pdf_mag[ic].dx), //mag
			      (pdf_magobs[ic].xmin+(float)(j)*pdf_magobs[ic].dx) -
			      (pdf_mag[ic].xmin+(float)(i)*pdf_mag[ic].dx), //mu
			      pdf_magobs[ic].y[j] * pdf_mag[ic].y[i]); //jpdf  
		    
		    /* Aparent distance modulus: mu_lambda = magobs - mag */
		    //Correcao de 2.0e-6 por causa da precisão do C
		    add2pdf(&pdf_mu[ic], (pdf_magobs[ic].xmin+(float)(j)*pdf_magobs[ic].dx) -
			    (pdf_mag[ic].xmin+(float)(i)*pdf_mag[ic].dx) + 2.0e-6, pdf_magobs[ic].y[j] * pdf_mag[ic].y[i]);
		  }
	      
	      for (imu=0; imu<pdf_mu[ic].n; imu++)
		//adoto o valor mu+(delta_mu/2) para evitar problema de precisão em add2pdf(pdf_mu0)
		mu[imu].fil[ic] = pdf_mu[ic].xmin + (float)(imu)*pdf_mu[ic].dx + (pdf_mu[ic].dx/2.0);
	      
	      if (print_pdfs==1) print_PDFmagfile (fp_mu, "PDFmu/", "pdf_mu_", name, "#fil mu PDF\n",
						   &pdf_mu[ic], 1.0e-6, ic);
	      nic++;
	      strcat(nomimag, nomi_mag[ic]) ;
	    } 
	}
      if (print_jpdfs==1) fclose (fp_jpdfmag);
      
      if (fit_av==1) // fitting extinction
	{
	  for (volume=0.0, iav=0; iav<pdf_av.n; iav++) // Av
	    {
	      /* Distance modulus in each filter: mu0_lambda = mu_lambda - AV(ALamb/AV) */
	      for (ic=1; ic<n_fil+1; ic++) // filtros
		{
		  set_pdf(&pdf_mu0_lamb[ic], 0, mu0_min, mu0_max, 0.005);
		  if (magobs.fil[ic] < -90.0 ||  smagobs.fil[ic] < -90.0) continue;
		  for (imu=0; imu<pdf_mu[ic].n; imu++) // mu0
		    add2pdf(&pdf_mu0_lamb[ic], mu[imu].fil[ic] -
			    ((pdf_av.xmin+(float)(iav)*pdf_av.dx)*alamb_av[ic]), pdf_mu[ic].y[imu]);
		}
	      
	      /* JPDF (mu0, Av): multiplies prob_mu0 of all filters for each Av */
	      for (imu=0; imu<pdf_mu0.n; imu++)
		{
		  jpdf_mu0_av[iav][imu]=1.0;
		  for (ic=1; ic<n_fil+1; ic++) 
		    {
		      if (magobs.fil[ic] < -90.0 || smagobs.fil[ic] < -90.0) continue;
		      jpdf_mu0_av[iav][imu] *= pdf_mu0_lamb[ic].y[imu];
		    }
		  if (magobs.fil[ic] < -90.0 || smagobs.fil[ic] < -90.0) continue;
		  
		  volume += jpdf_mu0_av[iav][imu];
		  pdf_av.y[iav] += jpdf_mu0_av[iav][imu]; //JPDF projection: PDF Av
		}
	    }      
	}
      else //use input extinction: NOT WORKING PROPERLY
	{
	  for (i=0; i<pdf_av.n; i++)
	    add2pdf(&pdf_av, -0.5+(float)(i)*pdf_av.dx+(avx+0.005),
		    gaussiana (1.0, savx, -0.5+(float)(i)*pdf_av.dx));
	  normalizes_pdf(&pdf_av);
	  
	  for (volume=0.0, iav=0; iav<pdf_av.n; iav++) // Av
	    {
	      /* Distance modulus in each filter: mu0_lambda = mu_lambda - AV(ALamb/AV) */
	      for (ic=1; ic<n_fil+1; ic++) // filtros
		{
		  set_pdf(&pdf_mu0_lamb[ic], 0, mu0_min, mu0_max, 0.005);
		  if (magobs.fil[ic] < -90.0 ||  smagobs.fil[ic] < -90.0) continue;
		  for (imu=0; imu<pdf_mu[ic].n; imu++) // mu0
		    add2pdf(&pdf_mu0_lamb[ic], mu[imu].fil[ic] -
			    ((pdf_av.xmin+(float)(iav)*pdf_av.dx)*alamb_av[ic]),
			    pdf_av.y[iav] * pdf_mu[ic].y[imu]);
		}
	      
	      /* JPDF (mu0, Av): multiplies prob_mu0 of all filters for each Av */
	      for (imu=0; imu<pdf_mu0.n; imu++)
		{
		  jpdf_mu0_av[iav][imu]=1.0;
		  for (ic=1; ic<n_fil+1; ic++) 
		    {
		      if (magobs.fil[ic] < -90.0 || smagobs.fil[ic] < -90.0) continue;
		      jpdf_mu0_av[iav][imu] *= pdf_mu0_lamb[ic].y[imu];
		    }
		  volume += jpdf_mu0_av[iav][imu];
		  pdf_av.y[iav] += jpdf_mu0_av[iav][imu]; //JPDF projection: PDF Av
		}
	    }
	}
      //printf("\nvolume %Lg", volume);

      if (nic != 0)
	{
	  /* PDF Av: median, mode and CI of 68%, 95% (JPDF projection) */ 
	  computeCI_pdf(&pdf_av) ;
	  print_PDFsummary (&pdf_av, "Av", "0");
	  if (print_pdfs==1) print_PDFfile (fp_av, "PDFav/", "pdf_av_", name, "#Av PDF\n",
					    &pdf_av, 1.0e-8, 0);
	  
	  /* PDF mu:median, mode and CI of 68%, 95% (JPDF projection) */
	  for (imu=0; imu<pdf_mu0.n; imu++)
	    for (iav=0; iav<pdf_av.n; iav++)
	      pdf_mu0.y[imu] += jpdf_mu0_av[iav][imu];
	  
	  computeCI_pdf(&pdf_mu0);
	  print_PDFsummary (&pdf_mu0, "mu0", "0");
	  if (print_pdfs==1) print_PDFfile (fp_mu0, "PDFmu0/", "pdf_mu0_", name, "#mu0 PDF\n",
					    &pdf_mu0, 1.0e-8, 0);
	  
	  /* Printing distance values */
	  if (pdf_mu0.ci_flag == 1)
	    {
	      printf("\n dist median = %.2f pc, CI(68%%): %.2f - %.2f, CI(95%%): %.2f - %.2f", 
		     10.0*pow(10.0,(0.2*pdf_mu0.median)),
		     10.0*pow(10.0,(0.2*pdf_mu0.me68min)), 10.0*pow(10.0,(0.2*pdf_mu0.me68max)),
		     10.0*pow(10.0,(0.2*pdf_mu0.me95min)), 10.0*pow(10.0,(0.2*pdf_mu0.me95max)));
	      printf("\n dist mode = %.2f pc, CI(68%%): %.2f - %.2f, CI(95%%): %.2f - %.2f\n", 
		     10.0*pow(10.0,(0.2*pdf_mu0.mode)),
		     10.0*pow(10.0,(0.2*pdf_mu0.mo68min)), 10.0*pow(10.0,(0.2*pdf_mu0.mo68max)),
		     10.0*pow(10.0,(0.2*pdf_mu0.mo95min)), 10.0*pow(10.0,(0.2*pdf_mu0.mo95max)));
	    }
	  else
	    {
	      printf("\n dist median = %.2f pc, CI(68%%): %.2f - %.2f, CI(95%%): %.2f - %.2f", 
		     ANNUL, ANNUL, ANNUL, ANNUL, ANNUL);
	      printf("\n dist mode = %.2f pc, CI(68%%): %.2f - %.2f, CI(95%%): %.2f - %.2f\n", 
		     ANNUL, ANNUL, ANNUL, ANNUL, ANNUL);	      
	    }

	  /* Normalizes JPDF(Av, mu0) */
	  float mu0_av_prob_min;
	  if (volume >= 1.0)
	    mu0_av_prob_min = 1.0e-6; 
	  else if (volume >= 1.0e-9 && volume < 1.0)
	    mu0_av_prob_min = 1.0e-7;
	  else if (volume >= 1.0e-19 && volume < 1.0e-9)
	    mu0_av_prob_min = 1.0e-10;
	  else mu0_av_prob_min = 1.0e-20;
	  for (ijpdf=0, iav=0; iav<pdf_av.n; iav++) 
	    for (imu=0; imu<pdf_mu0.n; imu++) 
	      //TSR: need to be improved
	      if ((jpdf_mu0_av[iav][imu]/volume) >= mu0_av_prob_min)
		{
		  /* Matrix (mu0, av, JPDF) */
		  mu0_av_prob[ijpdf][0] = pdf_mu0.xmin + (float)(imu)*pdf_mu0.dx ;
		  mu0_av_prob[ijpdf][1] = pdf_av.xmin + (float)(iav)*pdf_av.dx ;
		  mu0_av_prob[ijpdf][2] = jpdf_mu0_av[iav][imu]/volume;
		  ijpdf++;
	 	}  
	  //printf("\n IJPDF %d", ijpdf);
	  ijpdf--;
	  
	  /* Sorts matrix (mu0, Av, JPDF) in JPDF */
	  qsort( mu0_av_prob, ijpdf+1, sizeof mu0_av_prob[0], cmp );
	  
	  if (print_jpdfs==1)
	    {
	      strcpycat (temp, PDFSDIR, "JPDF_mu0_av/");
	      strcat (temp, "jpdf_mu0_av_");
	      strcat (temp, name);
	      fp_jpdfmu0av = fopen (temp, "w");
	      fprintf (fp_jpdfmu0av, "#mu0 Av JPDF\n");
	      for (i=0; i<ijpdf+1; i++)
		fprintf  (fp_jpdfmu0av, "%f %f %g\n", mu0_av_prob[i][0],
			  mu0_av_prob[i][1], mu0_av_prob[i][2]) ;
	      fclose (fp_jpdfmu0av) ;
	    }
	  
	  if (pdf_mu0.ci_flag == 1)
	    {
	      printf ("\n Mode JPDF(mu0,Av): %.2f, %.2f",
		      mu0_av_prob[ijpdf][0], mu0_av_prob[ijpdf][1]) ;
	      fprintf (fp_mo_out, "%.4f ", mu0_av_prob[ijpdf][1]);
	      print_outputfiles (fp_mo_out, fp_me_out, &pdf_av, 0) ;
	      fprintf (fp_mo_out, "%.4f ", mu0_av_prob[ijpdf][0]);
	      print_outputfiles (fp_mo_out, fp_me_out, &pdf_mu0, 0) ;
	      
	      fprintf (fp_me_out, "%.4f %.4f %.4f %.4f %.4f ",
		       10.0*pow(10.0,(0.2*pdf_mu0.median)),
		       10.0*pow(10.0,(0.2*pdf_mu0.me68min)), 10.0*pow(10.0,(0.2*pdf_mu0.me68max)),
		       10.0*pow(10.0,(0.2*pdf_mu0.me95min)), 10.0*pow(10.0,(0.2*pdf_mu0.me95max)));
	      fprintf (fp_mo_out, "%.4f %.4f %.4f %.4f %.4f %.4f ",
		       10.0*pow(10.0,(0.2*mu0_av_prob[ijpdf][0])),
		       10.0*pow(10.0,(0.2*pdf_mu0.mode)),
		       10.0*pow(10.0,(0.2*pdf_mu0.mo68min)), 10.0*pow(10.0,(0.2*pdf_mu0.mo68max)),
		       10.0*pow(10.0,(0.2*pdf_mu0.mo95min)), 10.0*pow(10.0,(0.2*pdf_mu0.mo95max)));
	    }
	  else
	    {
	      printf ("\n Mode JPDF(mu0,Av): %.2f, %.2f", ANNUL, ANNUL);
	      fprintf (fp_mo_out, "%.4f ", ANNUL);
	      print_outputfiles (fp_mo_out, fp_me_out, &pdf_av, 0) ;
	      fprintf (fp_mo_out, "%.4f ", ANNUL);
	      print_outputfiles (fp_mo_out, fp_me_out, &pdf_mu0, 0) ;
	      fprintf(fp_mo_out, "%.4f %.4f %.4f %.4f %.4f %.4f ",
		      ANNUL, ANNUL, ANNUL, ANNUL, ANNUL, ANNUL);
	      fprintf(fp_me_out, "%.4f %.4f %.4f %.4f %.4f ", ANNUL, ANNUL, ANNUL, ANNUL, ANNUL);
	    }
	  
	  /* Prints magobs used to calculate mu0 */
	  fprintf (fp_me_out, "%d %s\n", nic, nomimag) ;
	  fprintf (fp_mo_out, "%d %s\n", nic, nomimag) ;
	  printf("\n fil = %d: %s\n ", nic, nomimag) ;
	  
	  /* Distance PDF */
	  if (print_pdfs==1)
	    {
	      strcpycat (temp, PDFSDIR, "PDFdist/");
	      strcat (temp, "pdf_d_");
	      strcat (temp, name);
	      fp_dist = fopen (temp, "w");
	      fprintf (fp_dist, "#dist(pc) PDF\n" ); 
	      for (imu=0; imu<pdf_mu0.n; imu++)
		if (pdf_mu0.y[imu] >= 1.0e-8)
		  fprintf (fp_dist, "%f %g\n",
			   10.0*pow(10.0,(0.2*(pdf_mu0.xmin+(float)(imu)*pdf_mu0.dx))),
			   pdf_mu0.y[imu]*(5.0/log(10.0))/(10.0*pow(10.0,(0.2*(pdf_mu0.xmin+(float)(imu)*pdf_mu0.dx))))) ;
	      fclose (fp_dist) ;
	    }
	}
      else
	{
	  printf("\n Av median = %.2f, CI(68%%): %.2f - %.2f, CI(95%%): %.2f - %.2f", 
		 ANNUL, ANNUL, ANNUL, ANNUL, ANNUL);
	  printf("\n Av mode = %.2f, CI(68%%): %.2f - %.2f, CI(95%%): %.2f - %.2f", 
		 ANNUL, ANNUL, ANNUL, ANNUL, ANNUL);
	  printf("\n mu0 median = %.2f, CI(68%%): %.2f - %.2f, CI(95%%): %.2f - %.2f", 
		 ANNUL, ANNUL, ANNUL, ANNUL, ANNUL);
	  printf("\n mu0 mode = %.2f, CI(68%%): %.2f - %.2f, CI(95%%): %.2f - %.2f", 
		 ANNUL, ANNUL, ANNUL, ANNUL, ANNUL);
	  printf("\n dist median = %.2f pc, CI(68%%): %.2f - %.2f, CI(95%%): %.2f - %.2f", 
		 ANNUL, ANNUL, ANNUL, ANNUL, ANNUL);
	  printf("\n dist mode = %.2f pc, CI(68%%): %.2f - %.2f, CI(95%%): %.2f - %.2f\n", 
		 ANNUL, ANNUL, ANNUL, ANNUL, ANNUL);	  
	  printf ("\n Mode JPDF(mu0,Av): %.2f, %.2f", ANNUL, ANNUL);
	  
	  fprintf (fp_me_out, "%.4f %.4f %.4f %.4f %.4f ",
		   ANNUL, ANNUL, ANNUL, ANNUL, ANNUL);
	  fprintf (fp_mo_out, "%.4f %.4f %.4f %.4f %.4f %.4f ",
		   ANNUL, ANNUL, ANNUL, ANNUL, ANNUL, ANNUL);
	  
	  fprintf (fp_me_out, "%.4f %.4f %.4f %.4f %.4f ", 
		   ANNUL, ANNUL, ANNUL, ANNUL, ANNUL);
	  fprintf (fp_mo_out, "%.4f %.4f %.4f %.4f %.4f %.4f ", 
		   ANNUL, ANNUL, ANNUL, ANNUL, ANNUL, ANNUL);
	  
	  fprintf (fp_me_out, "%.4f %.4f %.4f %.4f %.4f ",
		   ANNUL, ANNUL, ANNUL, ANNUL, ANNUL);
	  fprintf (fp_mo_out, "%.4f %.4f %.4f %.4f %.4f %.4f ",
		   ANNUL, ANNUL, ANNUL, ANNUL, ANNUL, ANNUL);
	  
	  strcpy(nomimag,"null");
	  fprintf (fp_me_out, "%d %s\n", nic, nomimag) ;
	  fprintf (fp_mo_out, "%d %s\n", nic, nomimag) ;
	  printf("\n fil = %d: %s\n ", nic, nomimag) ;
	}
      
#ifdef ANDREAOUTPUT     
      // closes file with PDF
      fclose(fp_out_tmp);
#endif
      
      end = clock();
      printf("\nTime: %.2lf min\n",((double)(end - start)/CLOCKS_PER_SEC)/60.0);
  
    }
  
  fclose (fp_me_out) ;     
  fclose (fp_mo_out) ;
  if (print_input_pdfs==1)
    {
      fclose (fp_inppdfs_me); fclose (fp_inppdfs_mo);
    }
}
/************************************************************/


/* Coefficient extinction table */
void alamb_av_teff (float teff, float *alambda_av )
{
  int  i, ic;
  float fact;
  int nlin = 10, ncol = 17;
  float tab_alambda[10][17] = {
    //#Teff logg Av   Rv  Kepler    g      r      i     z  DDO51_finf  J     H      Ks    W1     W2     W3     W4
    {3500.0, 0.50, 1.0, 3.1, 0.718, 1.166, 0.861, 0.667, 0.485, 1.095, 0.292, 0.180, 0.119, 0.071, 0.055, 0.002, 0.000},
    {3750.0, 1.00, 1.0, 3.1, 0.757, 1.169, 0.863, 0.673, 0.487, 1.094, 0.293, 0.180, 0.118, 0.072, 0.056, 0.002, 0.000},
    {4000.0, 1.50, 1.0, 3.1, 0.780, 1.175, 0.865, 0.678, 0.488, 1.094, 0.293, 0.181, 0.119, 0.071, 0.055, 0.002, 0.000},
    {4250.0, 2.00, 1.0, 3.1, 0.797, 1.179, 0.866, 0.679, 0.489, 1.095, 0.292, 0.181, 0.119, 0.071, 0.056, 0.002, 0.000},
    {4500.0, 2.50, 1.0, 3.1, 0.809, 1.184, 0.867, 0.680, 0.489, 1.094, 0.293, 0.180, 0.119, 0.071, 0.055, 0.002, 0.000},
    {4750.0, 3.00, 1.0, 3.1, 0.818, 1.188, 0.868, 0.681, 0.491, 1.095, 0.294, 0.180, 0.118, 0.071, 0.056, 0.002, 0.000},
    {5000.0, 3.50, 1.0, 3.1, 0.828, 1.193, 0.868, 0.682, 0.490, 1.095, 0.294, 0.181, 0.118, 0.071, 0.055, 0.002, 0.000},
    {5250.0, 4.00, 1.0, 3.1, 0.836, 1.197, 0.869, 0.682, 0.491, 1.095, 0.294, 0.182, 0.118, 0.071, 0.055, 0.002, 0.000},
    {5500.0, 4.50, 1.0, 3.1, 0.843, 1.199, 0.870, 0.682, 0.491, 1.095, 0.294, 0.181, 0.118, 0.072, 0.055, 0.002, 0.000},
    {5750.0, 5.00, 1.0, 3.1, 0.850, 1.203, 0.870, 0.682, 0.492, 1.096, 0.295, 0.181, 0.118, 0.071, 0.056, 0.002, 0.000}
};

  if (teff<tab_alambda[0][0])
    i=0;
  else if (teff>=tab_alambda[nlin-1][0])
    i=nlin-2;
  else 
    for (i=0; i<nlin; i++)
      if (teff>=tab_alambda[i][0] && teff<tab_alambda[i+1][0])
	break;
  fact = (teff-tab_alambda[i][0])/(tab_alambda[i+1][0]-tab_alambda[i][0]);
  // for each filter:
  for (ic=4; ic<ncol; ic++)
    alambda_av[ic-3] = 
      (1.0-fact)*tab_alambda[i][ic] + fact*tab_alambda[i+1][ic] ;
}

/* Reading PARAM configuration */
void read_param_conf ()
{
  FILE *fp_conf;
  char	lixo[80];
  struct stat stt;
  
  if (stat ("param.conf", &stt) == 0)
    {
      printf ("\nReading param.conf...");
      fp_conf = fopen("param.conf", "r");
      fscanf (fp_conf, "%d", &metododireto);
      fgets(lixo, 80, fp_conf);
      fscanf (fp_conf, "%f %f %f %f %f %f", &use_teff, &use_lg, &use_dnu, &use_numax, &use_pspacing, &use_lum);
      fgets(lixo, 80, fp_conf);
      fscanf (fp_conf, "%d %d %d", &print_pdfs, &print_jpdfs, &print_input_pdfs);
      fgets(lixo, 80, fp_conf);
      fscanf(fp_conf, "%f", &max_smagobs);
      fgets(lixo, 80, fp_conf);
      fscanf(fp_conf, "%f", &delta_teff);
      fclose(fp_conf);
    }
  else
    {
      printf ("\nDefault param.conf:");
      metododireto = 0;
      use_teff = 1.0; use_lg = 0.0; use_dnu = 1.0; use_numax = 1.0; use_pspacing = 0.0; use_lum = 0.0;
      print_pdfs = 1; print_jpdfs = 0; print_input_pdfs = 0 ;
      max_smagobs = 0.1;
      delta_teff = 0.0;
    }
  (metododireto == 0) ? printf ("\n BAYESIAN METHOD") : printf ("\n DIRECT METHOD");
  printf("\n Using:");
  if (use_teff > 0) printf(" Teff");
  if (use_lg > 0) printf(" logg");
  if (use_dnu > 0) printf(" Dnu");
  if (use_numax > 0) printf(" numax");
  if (use_pspacing > 0) printf(" DPi1");
  if (use_lum > 0) printf(" L");
  if (print_pdfs==1) printf("\n Will print PDFs...");
  if (print_jpdfs==1) printf("\n Will print JPDFs...");
  if (print_input_pdfs==1) printf("\n Will print input PDFs...");
  printf ("\n Maximum error in apparent magnitudes to compute distances/extinctions = %.2f mag", max_smagobs);
  if (delta_teff > 0.0) printf ("\n Will apply an offset in Teff: Delta_teff = %.2f K", delta_teff); 
  printf("\n");
  if ((print_pdfs==1) || (print_jpdfs==1) || (print_input_pdfs==1)) makedir (PDFSDIR);
}
