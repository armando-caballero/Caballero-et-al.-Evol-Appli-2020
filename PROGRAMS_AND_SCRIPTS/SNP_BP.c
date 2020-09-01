// SNP_BP_DRUET.c

#include "libhdr"

#define NN 50000 // Maximum 50000 SNPs segregating
#define CC 20001   // Maximum N=10000 individuals
#define MC 20000 // Maximum 20000 SNPs per chromosome

int i, j, s, m, x, nind, HOM1, HOM2, trait, crossover, rnd;
int mut[NN], pos[NN], crom[CC][MC], num_mut_crom[CC];
int numSNP, numQTL;
double w, ps[NN], a[NN], h[NN], q[NN], FITN[CC], FITN_F[CC], PHEN[CC], PHEN_F[CC], HML[CC], HMLmaf[CC], HMLqtl[CC], VE;
double randE[CC], scale_a, MAF;
double IDFhatI, IDFhatII, IDFhatIII, IDFhom, IDFexh, IDFroh1, IDFroh100, IDFroh1000;
double SqE_FhatI, SqE_FhatII, SqE_FhatIII, SqE_Fhom, SqE_Fexh, SqE_Froh1, SqE_Froh100, SqE_Froh1000;
double d_a, ID_2dpq, sum_aq, mean_FITN, mean_FITN_F, mean_PHEN, mean_PHEN_F;
char ch;

struct acc Phe, AHML, AHMLmaf, AHMLqtl, AFhatI, AFhatII, AFhatIII, AFhom, AFexh, AFroh1, AFroh100, AFroh1000;
struct acc Ave_q_qtl, Ave_q_snp;

struct covacc PheHML, PheHMLmaf, PheHMLqtl, PheFhatI, PheFhatII, PheFhatIII, PheFhom, PheFexh, PheFroh1, PheFroh100, PheFroh1000;
struct covacc HMLHMLmaf, HMLHMLqtl, HMLFhatI, HMLFhatII, HMLFhatIII, HMLFhom, HMLFexh, HMLFroh1, HMLFroh100, HMLFroh1000;
struct covacc HMLmafHMLqtl, HMLmafFhatI, HMLmafFhatII, HMLmafFhatIII, HMLmafFhom, HMLmafFexh, HMLmafFroh1, HMLmafFroh100, HMLmafFroh1000;
struct covacc HMLqtlFhatI, HMLqtlFhatII, HMLqtlFhatIII, HMLqtlFhom, HMLqtlFexh, HMLqtlFroh1, HMLqtlFroh100, HMLqtlFroh1000;
struct covacc FhatIFhatII, FhatIFhatIII, FhatIFhom, FhatIFexh, FhatIFroh1, FhatIFroh100, FhatIFroh1000;
struct covacc FhatIIFhatIII, FhatIIFhom, FhatIIFexh, FhatIIFroh1, FhatIIFroh100, FhatIIFroh1000;
struct covacc FhatIIIFhom, FhatIIIFexh, FhatIIIFroh1, FhatIIIFroh100, FhatIIIFroh1000;
struct covacc FhomFexh, FhomFroh1, FhomFroh100, FhomFroh1000;
struct covacc FexhFroh1, FexhFroh100, FexhFroh1000;
struct covacc Froh1Froh100, Froh1Froh1000;
struct covacc Froh100Froh1000;

struct acc q_SNP[201], q_QTL[201], a_QTL[201], ID_QTL[201];

FILE *fdat, *fped, *fout, *fmap, *fPphen, *fallsnps, *fLISTQTL;
FILE *fFfile, *foutfile, *fsumID, *fsumM, *fsumV, *fsumR, *fsumSqE, *fcheck;
FILE *fdis;

main()
{
	tracestart();
	getseed();
	getintandskip("Trait (0; fitness, 1: QT) :",&trait,0,1);
	getintandskip("NIND :",&nind,10,10000);
	getrealandskip("scale_a :",&scale_a,0.0,10.0);
	getrealandskip("VE :",&VE,0.0,100.0);
	getrealandskip("MAF for HML :",&MAF,0.0,1.0);

	readfiles();
	distribution();
	PLINK_files();
	int status = system("bash shell_F_values");
	regression_pm_F();
	printing();
	writeseed();

	return(0);
}

/* **************************************************************************** */

readfiles()
{
	fout = fopen ("outfileSLIM","w");
	fsumID = fopen ("summaryID","w");
//	fsumM = fopen ("summaryM","w");
//	fsumV = fopen ("summaryV","w");
	fsumR = fopen ("summaryR","w");
	fsumSqE = fopen ("summarySqE","w");
	foutfile = fopen ("outfile","w");
	fcheck = fopen ("checkfile","w");
	fdis = fopen ("distributionfile","w");


	fallsnps = fopen ("list_allsnps","w");
	fLISTQTL = fopen ("list_qtls","w");

	// ********** Read slimout to get SNP positions, effects and frequencies ********** 

	fdat = fopen ("slimout","r");

	// ********** gets the position, effects and frequencies of all mutations simulated (numSNP) **********

	lookfortext("Mutations:");

	while (!feof(fdat))
	{
		s ++;
		fscanf(fdat,"%d", &x);
		mut[s] = x;
		for (j=1; j<=4; j++)	fscanf(fdat,"%c", &ch);
		fscanf(fdat,"%d", &x);
		pos[s] = x;
		fscanf(fdat,"%lf", &w);
//		if (w < -1.0)	w=-1.0;
		ps[s] = w;
		fscanf(fdat,"%lf", &w);
		if (w == -99.0) w=0.0;
		if (trait == 0) a[s] = ps[s];
		else
		{
			if (w < 0.0)	a[s] = w * scale_a;
			else			a[s] = -w * scale_a;
		}
		if (a[s] != 0.0)	numQTL ++;
		fscanf(fdat,"%lf", &w);
		h[s] = w;
		if (ps[s]==0.0) h[s] = 0.0;
		fscanf(fdat,"%d", &x);
		q[s] = x / (2.0*nind);

		if (a[s] != 0.0)	accum (&Ave_q_qtl, q[s]);
		else				accum (&Ave_q_snp, q[s]);

		fscanf(fdat,"%c", &ch);
		fscanf(fdat,"%c", &ch);
		if (ch != 'G')	ungetc(ch, fdat);
		else		break;

		// Inbreeding Depression
		if (a[s] != 0.0)
		{
			d_a = a[s] * (h[s] - 0.5);
			sum_aq += a[s] * q[s];
			ID_2dpq += (2.0 * d_a * q[s] * (1.0 - q[s]));
		}
	}
	numSNP = s;

	fclose(fdat);

	// ********** gets mutations **********

	fdat = fopen ("slimout","r");
	lookfortext("Genomes:");

	fscanf(fdat,"%c", &ch);

	for (i=1; i<=(2*nind);i++)
	{
		m = 0;

		fscanf(fdat,"%c", &ch);
		fscanf(fdat,"%c", &ch);
		fscanf(fdat,"%c", &ch);
		fscanf(fdat,"%d", &x);

		while (!feof(fdat))
		{
			fscanf(fdat,"%c", &ch);

			if (ch == '\n')	break;
			else
			{
				ungetc(ch, fdat);
				m ++;
				fscanf(fdat,"%d", &x);
				crom[i][m] = x;
			}			
		}
		num_mut_crom[i] = m;
		if (tracelevel!=0)  if (m > MC)   fprintf(fcheck, "*****num_mut_crom > MAX\n");
		if (m > MC)   fprintf(foutfile, "*****num_mut_crom > MAX\n");
	}

	fclose(fdat);

	return(0);
}

/* **************************************************************************** */

distribution()
{
	for (s=1; s<=numSNP; s++)
	{
		if (a[s] != 0.0) 
		{
			for (i=0; i<=100; i++)
			if ( (q[s] >= i/100.0) && (q[s] < (i+1)/100.0) )
			{
				accum (&q_QTL[i], 1.0);
				accum (&a_QTL[i], a[s]);
				d_a = a[s] * (h[s] - 0.5);
				accum (&ID_QTL[i], 2.0 * d_a * q[s] * (1.0 - q[s]));
			}
		}
		else
		{
			for (i=0; i<=100; i++)
			if ( (q[s] >= i/100.0) && (q[s] < (i+1)/100.0) )
			accum (&q_SNP[i], 1.0);
		}
	}

	for (i=0; i<=100; i++)
	fprintf (fdis, "%4.2f %6.2f %6.2f %f %f\n", i/100.0, accsum(&q_SNP[i]), accsum(&q_QTL[i]), accmean(&a_QTL[i]), accsum(&ID_QTL[i]));
}

/* **************************************************************************** */

PLINK_files()
{
	fped = fopen ("data.ped","w");
	fmap = fopen ("data.map","w");
	fPphen = fopen ("qt.phe","w");

	fped = fopen ("data.ped","a");
	fmap = fopen ("data.map","a");
	fPphen = fopen ("qt.phe","a");

	double sum, sum2;

	for (i=1; i<(2*nind);i+=2)
	{
		if (i%2 != 0)
		{
			fprintf(fped,"1 IND%d 0 0 1 -9 ", (i+1)/2);
			FITN[(i+1)/2] = 1.0;
			FITN_F[(i+1)/2] = 1.0;
			PHEN[(i+1)/2] = 0.0;
			PHEN_F[(i+1)/2] = 0.0;
		}

		if (uniform() < 0.5) rnd = 1;
		else				 rnd = 2;

		if (tracelevel!=0)     fprintf(fcheck, "\nIND %d\n", (i+1)/2);

		for (s=1; s<=numSNP; s++)
		{
			if (i == 1)
			{
				// PLINK POSFILE
				if (a[s] == 0.0)	fprintf(fmap,"1 SNP%d 0 %d\n", s, pos[s]);
			}
			HOM1 = 0; HOM2 = 0;

			for (m=1; m<=num_mut_crom[i]; m++)
			{
				if (crom[i][m] == mut[s])
				{
					HOM1 = 1;
					break;
				}
			}
			for (m=1; m<=num_mut_crom[i+1]; m++)
			{
				if (crom[i+1][m] == mut[s])
				{
					HOM2 = 1;
					break;
				}
			}

			if (tracelevel!=0) if (a[s] != 0.0) fprintf(fcheck, "s = %d  HOM1 = %d  HOM2 = %d  ", s, HOM1, HOM2);

			if ((HOM1==0) && (HOM2==0))
			{
				if (a[s] == 0.0)
				{
					fprintf(fped,"A A ");
//					if ((1.0-q[s]) <= MAF)	HMLmaf[(i+1)/2]++;
//					HML[(i+1)/2]++;
				}
			}
			else if ((HOM1==1) && (HOM2==0))
			{
				if (a[s] == 0.0)	fprintf(fped,"T A ");
				if (trait == 0)		FITN[(i+1)/2] *= (1.0 + ps[s]*h[s]);
				else if (a[s] != -99)	PHEN[(i+1)/2] += a[s]*h[s];

				if (rnd == 1)
				{
					if (trait == 0)		FITN_F[(i+1)/2] *= (1.0 + ps[s]);
					else if (a[s] != -99)	PHEN_F[(i+1)/2] += a[s];
				}
				else		/*11*/;
			}
			else if ((HOM1==0) && (HOM2==1))
			{
				if (a[s] == 0.0)	fprintf(fped,"A T ");
				if (trait == 0)		FITN[(i+1)/2] *= (1.0 + ps[s]*h[s]);	
				else if (a[s] != 0)		PHEN[(i+1)/2] += a[s]*h[s];

				if (rnd == 1)	/*11*/;
				else
				{
					if (trait == 0)		FITN_F[(i+1)/2] *= (1.0 + ps[s]);	
					else if (a[s] != 0)		PHEN_F[(i+1)/2] += a[s];
				}
			}
			else
			{
				if (a[s] == 0.0)
				{
					fprintf(fped,"T T ");
					if (q[s] <= MAF)	HMLmaf[(i+1)/2]++;
					HML[(i+1)/2]++;
				}
				if (trait == 0)		FITN[(i+1)/2] *= (1.0 + ps[s]);
				else if (a[s] != 0)		PHEN[(i+1)/2] += a[s];

				if (trait == 0)		FITN_F[(i+1)/2] *= (1.0 + ps[s]);
				else if (a[s] != 0)		PHEN_F[(i+1)/2] += a[s];

				if (a[s] != 0.0)	HMLqtl[(i+1)/2]++;
			}

			if (tracelevel!=0) if (a[s] != 0.0) fprintf(fcheck, "PHEN = %f    PHEN_F = %f\n", PHEN[(i+1)/2], PHEN_F[(i+1)/2]);
		}
		fprintf(fped,"\n");

		if (i%2 != 0)
		{
			randE[(i+1)/2] = normal(0.0,VE);
			// PLINK PHENOFILE
			if (trait == 0)	fprintf(fPphen,"%f %f %f %f\n", log(FITN[(i+1)/2]), HML[(i+1)/2], HMLmaf[(i+1)/2], HMLqtl[(i+1)/2]);
			else				fprintf(fPphen,"%f %f %f %f\n", PHEN[(i+1)/2] + randE[(i+1)/2], HML[(i+1)/2], HMLmaf[(i+1)/2], HMLqtl[(i+1)/2]);

			// VARIANCE OF GENOTYPIC VALUES FOR THE TRAIT
			sum += PHEN[(i+1)/2];
			sum2 += (PHEN[(i+1)/2] * PHEN[(i+1)/2]);
		}
	}

	// List all SNPs
//	fprintf(fallsnps,"SNP  pos  s  a  h  q\n");
	fprintf(fallsnps,"%d\n", numSNP);
	for (s=1; s<=numSNP; s++)	fprintf(fallsnps,"%d  %d  %f  %f  %f  %f\n", s, pos[s], ps[s], a[s], h[s], q[s]);

	//FILE WITH QTLs
	fprintf(fLISTQTL,"SNP  pos  s  a  h  q\n", s, pos[s], ps[s], a[s], h[s], q[s]);
	for (s=1; s<=numSNP; s++)
	if (a[s] != 0)
	fprintf(fLISTQTL,"%d  %d  %f  %f  %f  %f\n", s, pos[s], ps[s], a[s], h[s], q[s]);

	fclose(fped);
	fclose(fmap);
	fclose(fPphen);

	return(0);
}

/* **************************************************************************** */

regression_pm_F ()
{
	double w, genvalue[CC], HomoMutLoad[CC], HomoMutLoadmaf[CC], HomoMutLoadqtl[CC], FhatI[CC], FhatII[CC], FhatIII[CC], Fhom[CC], Fexh[CC], Froh1[CC], Froh100[CC], Froh1000[CC];

    fFfile = fopen ("data.F","r");

	for (i=1; i<=nind; i++)
	{
		fscanf(fFfile,"%lf", &w);
		genvalue[i] = w;

		fscanf(fFfile,"%lf", &w);
		HomoMutLoad[i] = w;

		fscanf(fFfile,"%lf", &w);
		HomoMutLoadmaf[i] = w;

		fscanf(fFfile,"%lf", &w);
		HomoMutLoadqtl[i] = w;

		fscanf(fFfile,"%lf", &w);
		FhatI[i] = w;

		fscanf(fFfile,"%lf", &w);
		FhatII[i] = w;

		fscanf(fFfile,"%lf", &w);
		FhatIII[i] = w;

		fscanf(fFfile,"%lf", &w);
		Fhom[i] = w;

		fscanf(fFfile,"%lf", &w);
		Fexh[i] = w;

		fscanf(fFfile,"%lf", &w);
		Froh1[i] = w;

		fscanf(fFfile,"%lf", &w);
		Froh100[i] = w;

		fscanf(fFfile,"%lf", &w);
		Froh1000[i] = w;
	}

	close(fFfile);

	for (i=1; i<=nind; i++)
	{
		accum (&Phe, genvalue[i]);
		accum (&AHML, HomoMutLoad[i]);
		accum (&AHMLmaf, HomoMutLoadmaf[i]);
		accum (&AHMLqtl, HomoMutLoadqtl[i]);
		accum (&AFhatI, FhatI[i]);
		accum (&AFhatII, FhatII[i]);
		accum (&AFhatIII, FhatIII[i]);
		accum (&AFhom, Fhom[i]);
		accum (&AFexh, Fexh[i]);
		accum (&AFroh1, Froh1[i]);
		accum (&AFroh100, Froh100[i]);
		accum (&AFroh1000, Froh1000[i]);

		covaccum (&PheHML, genvalue[i], HomoMutLoad[i]);
		covaccum (&PheHMLmaf, genvalue[i], HomoMutLoadmaf[i]);
		covaccum (&PheHMLqtl, genvalue[i], HomoMutLoadqtl[i]);
		covaccum (&PheFhatI, genvalue[i], FhatI[i]);
		covaccum (&PheFhatII, genvalue[i], FhatII[i]);
		covaccum (&PheFhatIII, genvalue[i], FhatIII[i]);
		covaccum (&PheFhom, genvalue[i], Fhom[i]);
		covaccum (&PheFexh, genvalue[i], Fexh[i]);
		covaccum (&PheFroh1, genvalue[i], Froh1[i]);
		covaccum (&PheFroh100, genvalue[i], Froh100[i]);
		covaccum (&PheFroh1000, genvalue[i], Froh1000[i]);

		covaccum (&HMLHMLmaf, HomoMutLoad[i], HomoMutLoadmaf[i]);
		covaccum (&HMLHMLqtl, HomoMutLoad[i], HomoMutLoadqtl[i]);
		covaccum (&HMLFhatI, HomoMutLoad[i], FhatI[i]);
		covaccum (&HMLFhatII, HomoMutLoad[i], FhatII[i]);
		covaccum (&HMLFhatIII, HomoMutLoad[i], FhatIII[i]);
		covaccum (&HMLFhom, HomoMutLoad[i], Fhom[i]);
		covaccum (&HMLFexh, HomoMutLoad[i], Fexh[i]);
		covaccum (&HMLFroh1, HomoMutLoad[i], Froh1[i]);
		covaccum (&HMLFroh100, HomoMutLoad[i], Froh100[i]);
		covaccum (&HMLFroh1000, HomoMutLoad[i], Froh1000[i]);

		covaccum (&HMLmafHMLqtl, HomoMutLoadmaf[i], HomoMutLoadqtl[i]);
		covaccum (&HMLmafFhatI, HomoMutLoadmaf[i], FhatI[i]);
		covaccum (&HMLmafFhatII, HomoMutLoadmaf[i], FhatII[i]);
		covaccum (&HMLmafFhatIII, HomoMutLoadmaf[i], FhatIII[i]);
		covaccum (&HMLmafFhom, HomoMutLoadmaf[i], Fhom[i]);
		covaccum (&HMLmafFexh, HomoMutLoadmaf[i], Fexh[i]);
		covaccum (&HMLmafFroh1, HomoMutLoadmaf[i], Froh1[i]);
		covaccum (&HMLmafFroh100, HomoMutLoadmaf[i], Froh100[i]);
		covaccum (&HMLmafFroh1000, HomoMutLoadmaf[i], Froh1000[i]);

		covaccum (&HMLqtlFhatI, HomoMutLoadqtl[i], FhatI[i]);
		covaccum (&HMLqtlFhatII, HomoMutLoadqtl[i], FhatII[i]);
		covaccum (&HMLqtlFhatIII, HomoMutLoadqtl[i], FhatIII[i]);
		covaccum (&HMLqtlFhom, HomoMutLoadqtl[i], Fhom[i]);
		covaccum (&HMLqtlFexh, HomoMutLoadqtl[i], Fexh[i]);
		covaccum (&HMLqtlFroh1, HomoMutLoadqtl[i], Froh1[i]);
		covaccum (&HMLqtlFroh100, HomoMutLoadqtl[i], Froh100[i]);
		covaccum (&HMLqtlFroh1000, HomoMutLoadqtl[i], Froh1000[i]);

		covaccum (&FhatIFhatII, FhatI[i], FhatII[i]);
		covaccum (&FhatIFhatIII, FhatI[i], FhatIII[i]);
		covaccum (&FhatIFhom, FhatI[i], Fhom[i]);
		covaccum (&FhatIFexh, FhatI[i], Fexh[i]);
		covaccum (&FhatIFroh1, FhatI[i], Froh1[i]);
		covaccum (&FhatIFroh100, FhatI[i], Froh100[i]);
		covaccum (&FhatIFroh1000, FhatI[i], Froh1000[i]);

		covaccum (&FhatIIFhatIII, FhatII[i], FhatIII[i]);
		covaccum (&FhatIIFhom, FhatII[i], Fhom[i]);
		covaccum (&FhatIIFexh, FhatII[i], Fexh[i]);
		covaccum (&FhatIIFroh1, FhatII[i], Froh1[i]);
		covaccum (&FhatIIFroh100, FhatII[i], Froh100[i]);
		covaccum (&FhatIIFroh1000, FhatII[i], Froh1000[i]);

		covaccum (&FhatIIIFhom, FhatIII[i], Fhom[i]);
		covaccum (&FhatIIIFexh, FhatIII[i], Fexh[i]);
		covaccum (&FhatIIIFroh1, FhatIII[i], Froh1[i]);
		covaccum (&FhatIIIFroh100, FhatIII[i], Froh100[i]);
		covaccum (&FhatIIIFroh1000, FhatIII[i], Froh1000[i]);

		covaccum (&FhomFexh, Fhom[i], Fexh[i]);
		covaccum (&FhomFroh1, Fhom[i], Froh1[i]);
		covaccum (&FhomFroh100, Fhom[i], Froh100[i]);
		covaccum (&FhomFroh1000, Fhom[i], Froh1000[i]);

		covaccum (&FexhFroh1, Fexh[i], Froh1[i]);
		covaccum (&FexhFroh100, Fexh[i], Froh100[i]);
		covaccum (&FexhFroh1000, Fexh[i], Froh1000[i]);

		covaccum (&Froh1Froh100, Froh1[i], Froh100[i]);
		covaccum (&Froh1Froh1000, Froh1[i], Froh1000[i]);

		covaccum (&Froh100Froh1000, Froh100[i], Froh1000[i]);
	}

	IDFhatI = covariance(&PheFhatI) / variance(&AFhatI);
	IDFhatII = covariance(&PheFhatII) / variance(&AFhatII);
	IDFhatIII = covariance(&PheFhatIII) / variance(&AFhatIII);
	IDFhom = covariance(&PheFhom) / variance(&AFhom);
	IDFexh = covariance(&PheFexh) / variance(&AFexh);
	IDFroh1 = covariance(&PheFroh1) / variance(&AFroh1);
	IDFroh100 = covariance(&PheFroh100) / variance(&AFroh100);
	IDFroh1000 = covariance(&PheFroh1000) / variance(&AFroh1000);

	SqE_FhatI = pow(ID_2dpq + IDFhatI, 2.0);
	SqE_FhatII = pow(ID_2dpq + IDFhatII, 2.0);
	SqE_FhatIII = pow(ID_2dpq + IDFhatIII, 2.0);
	SqE_Fhom = pow(ID_2dpq + IDFhom, 2.0);
	SqE_Fexh = pow(ID_2dpq + IDFexh, 2.0);
	SqE_Froh1 = pow(ID_2dpq + IDFroh1, 2.0);
	SqE_Froh100 = pow(ID_2dpq + IDFroh100, 2.0);
	SqE_Froh1000 = pow(ID_2dpq + IDFroh1000, 2.0);
}


/* ***************************************************************** */


printing()
{
	int hom;

	fprintf(fout,"Positions of SNPs\n\n");
	for (s=1; s<=numSNP; s++)
	{
		fprintf(fout,"SNP%d  mut=%d  pos=%d  ps=%f  a=%f  h=%f  q=%f\n", s, mut[s], pos[s], ps[s], a[s], h[s], q[s]);
	}
	fprintf(fout,"\n");

	for (i=1; i<=(2*nind);i++)
	{
		fprintf(fout,"CHROM %d  ", i);

		for (m=1; m<=num_mut_crom[i]; m++)
		fprintf(fout,"%d ", crom[i][m]);

		fprintf(fout,"\n");
	}

/*	fprintf(fout,"\n\Chromosome of parents\n\n");
	for (i=1; i<=nind;i++)
	for (hom=1; hom<=2;hom++)
	{
		fprintf(fout,"ind %d chromosomeP %d   ", i, hom);
		for (s=1; s<=numSNP; s++)
		{
			fprintf(fout,"%d ", chromosomeP[i][s][hom]);
		}
		fprintf(fout,"\n");
	}
*/
	fprintf(fout,"\n\Phenotypes of parents and parents_F=1\n\n");
	for (i=1; i<=nind;i++)	fprintf(fout,"%d %f %f\n", i, PHEN[i], PHEN_F[i]);

//	outfile

	for (i=1; i<nind; i++)
	{
		mean_FITN += FITN[i]/(double)nind;
		mean_FITN_F += FITN_F[i]/(double)nind;
		mean_PHEN += PHEN[i]/(double)nind;
		mean_PHEN_F += PHEN_F[i]/(double)nind;
	}

	if (tracelevel!=0)     fprintf(fcheck, "\nmeanPHEN=%f   meanPHEN_F=%f\n", mean_PHEN, mean_PHEN_F);


	fprintf(foutfile,"\nN=%d  NSNPs=%d  NQTLs=%d  MAF_HML=%f  Ave_q_qtl=%f sd=%f  Ave_q_snp=%f sd=%f\n", nind, numSNP, numQTL, MAF, accmean(&Ave_q_qtl), se(&Ave_q_qtl), sqrt(variance(&Ave_q_snp)), sqrt(variance(&Ave_q_snp)));
	fprintf(foutfile, "AFhatI=%6.4f    AFhatII=%6.4f    AFhatIII=%6.4f\n", accmean(&AFhatI), accmean(&AFhatII), accmean(&AFhatIII)); 
	fprintf(foutfile, "AFhom=%6.4f    AFexh=%6.4f\n", accmean(&AFhom), accmean(&AFexh)); 
	fprintf(foutfile, "AFroh1=%6.4f    AFroh100=%6.4f    AFroh1000=%6.4f\n", accmean(&AFroh1), accmean(&AFroh100), accmean(&AFroh1000)); 

	fprintf(foutfile, "VFhatI=%6.4f    VFhatII=%6.4f    VFhatIII=%6.4f\n", variance(&AFhatI), variance(&AFhatII), variance(&AFhatIII)); 
	fprintf(foutfile, "VFhom=%6.4f    VFexh=%6.4f\n", variance(&AFhom), variance(&AFexh)); 
	fprintf(foutfile, "VFroh1=%6.4f    VFroh100=%6.4f    VFroh1000=%6.4f\n", variance(&AFroh1), variance(&AFroh100), variance(&AFroh1000)); 

	if (trait == 0)	fprintf(foutfile, "mean_FITN=%6.4f    mean_FITN_F=%6.4f    ID_F1=%6.4f\n", mean_FITN, mean_FITN_F, -log(mean_FITN_F/mean_FITN)); 
	else			fprintf(foutfile, "mean_PHEN=%6.4f    ExpPHEN=%6.4f    mean_PHEN_F=%6.4f    ID_F1=%6.4f\n", mean_PHEN, sum_aq+ID_2dpq, mean_PHEN_F, mean_PHEN-mean_PHEN_F); 

	fprintf(foutfile, "            HML      HMLmaf   HMLqtl   FhatI    FhatII   FhatIII   Fhom    Fexh     Froh1   Froh100   Froh1000\n");
	fprintf(foutfile, "Phe        %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
		correlation(&PheHML), correlation(&PheHMLmaf), correlation(&PheHMLqtl), correlation(&PheFhatI), correlation(&PheFhatII), correlation(&PheFhatIII), correlation(&PheFhom), correlation(&PheFexh), correlation(&PheFroh1), correlation(&PheFroh100), correlation(&PheFroh1000)); 
	fprintf(foutfile, "HML                  %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f\n",
		correlation(&HMLHMLmaf), correlation(&HMLHMLqtl), correlation(&HMLFhatI), correlation(&HMLFhatII), correlation(&HMLFhatIII), correlation(&HMLFhom), correlation(&HMLFexh), correlation(&HMLFroh1), correlation(&HMLFroh100), correlation(&HMLFroh1000)); 
	fprintf(foutfile, "HMLmaf                        %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f\n",
		correlation(&HMLmafHMLqtl), correlation(&HMLmafFhatI), correlation(&HMLmafFhatII), correlation(&HMLmafFhatIII), correlation(&HMLmafFhom), correlation(&HMLmafFexh), correlation(&HMLmafFroh1), correlation(&HMLmafFroh100), correlation(&HMLmafFroh1000)); 
	fprintf(foutfile, "HMLqtl                                 %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f\n",
		correlation(&HMLqtlFhatI), correlation(&HMLqtlFhatII), correlation(&HMLqtlFhatIII), correlation(&HMLqtlFhom), correlation(&HMLqtlFexh), correlation(&HMLqtlFroh1), correlation(&HMLqtlFroh100), correlation(&HMLqtlFroh1000)); 

	fprintf(foutfile, "FhatI                                           %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f\n",
		correlation(&FhatIFhatII), correlation(&FhatIFhatIII), correlation(&FhatIFhom), correlation(&FhatIFexh), correlation(&FhatIFroh1), correlation(&FhatIFroh100), correlation(&FhatIFroh1000)); 
	fprintf(foutfile, "FhatII                                                   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f\n",
		correlation(&FhatIIFhatIII), correlation(&FhatIIFhom), correlation(&FhatIIFexh), correlation(&FhatIIFroh1), correlation(&FhatIIFroh100), correlation(&FhatIIFroh1000)); 
	fprintf(foutfile, "FhatIII                                                           %6.4f   %6.4f   %6.4f   %6.4f   %6.4f\n",
		correlation(&FhatIIIFhom), correlation(&FhatIIIFexh), correlation(&FhatIIIFroh1), correlation(&FhatIIIFroh100), correlation(&FhatIIIFroh1000)); 
	fprintf(foutfile, "Fhom                                                                       %6.4f   %6.4f   %6.4f   %6.4f\n",
		correlation(&FhomFexh), correlation(&FhomFroh1), correlation(&FhomFroh100), correlation(&FhomFroh1000)); 
	fprintf(foutfile, "Fexh                                                                                %6.4f   %6.4f   %6.4f\n",
		correlation(&FexhFroh1), correlation(&FexhFroh100), correlation(&FexhFroh1000)); 
	fprintf(foutfile, "Froh1                                                                                        %6.4f   %6.4f\n", correlation(&Froh1Froh100), correlation(&Froh1Froh1000)); 
	fprintf(foutfile, "Froh100                                                                                               %6.4f\n", correlation(&Froh100Froh1000));

	fprintf(foutfile, "ID_2dpq=%6.4f\n", ID_2dpq); 
	fprintf(foutfile, "IDFhatI=%6.4f    IDFhatII=%6.4f    IDFhatIII=%6.4f\n", -IDFhatI, -IDFhatII, -IDFhatIII); 
	fprintf(foutfile, "IDFhom=%6.4f    IDFexh=%6.4f\n", -IDFhom, -IDFexh); 
	fprintf(foutfile, "IDFroh1=%6.4f    IDFroh100=%6.4f    IDFroh1000=%6.4f\n", -IDFroh1, -IDFroh100, -IDFroh1000); 

//	Summary ID file

	if (trait == 0)	fprintf(fsumID, "%6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n", ID_2dpq, mean_FITN-mean_FITN_F, -IDFhatI, -IDFhatII, -IDFhatIII, -IDFhom, -IDFexh, -IDFroh1, -IDFroh100, -IDFroh1000); 
	else			fprintf(fsumID, "%6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n", ID_2dpq, mean_PHEN-mean_PHEN_F, -IDFhatI, -IDFhatII, -IDFhatIII, -IDFhom, -IDFexh, -IDFroh1, -IDFroh100, -IDFroh1000); 

//	Summary Mean file

//	fprintf(fsumM, "%6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n", accmean(&AFhatI), accmean(&AFhatII), accmean(&AFhatIII), accmean(&AFhom), accmean(&AFexh), accmean(&AFroh1), accmean(&AFroh100), accmean(&AFroh1000)); 

//	Summary Var file

//	fprintf(fsumV, "%6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n", variance(&AFhatI), variance(&AFhatII), variance(&AFhatIII), variance(&AFhom), variance(&AFexh), variance(&AFroh1), variance(&AFroh100), variance(&AFroh1000)); 

//	Summary SqE file

	if (trait == 0)	fprintf(fsumSqE, "%6.4f  %6.4f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f\n", ID_2dpq, mean_FITN-mean_FITN_F, SqE_FhatI, SqE_FhatII, SqE_FhatIII, SqE_Fhom, SqE_Fexh, SqE_Froh1, SqE_Froh100, SqE_Froh1000); 
	else				fprintf(fsumSqE, "%6.4f  %6.4f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f\n", ID_2dpq, mean_PHEN-mean_PHEN_F, SqE_FhatI, SqE_FhatII, SqE_FhatIII, SqE_Fhom, SqE_Fexh, SqE_Froh1, SqE_Froh100, SqE_Froh1000); 

//	SummaryR (correlations) file

	fprintf(fsumR, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f ",
		correlation(&PheHML), correlation(&PheHMLmaf), correlation(&PheHMLqtl), correlation(&PheFhatI), correlation(&PheFhatII), correlation(&PheFhatIII), correlation(&PheFhom), correlation(&PheFexh), correlation(&PheFroh1), correlation(&PheFroh100), correlation(&PheFroh1000)); 
	fprintf(fsumR, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f ",
		correlation(&HMLHMLmaf), correlation(&HMLHMLqtl), correlation(&HMLFhatI), correlation(&HMLFhatII), correlation(&HMLFhatIII), correlation(&HMLFhom), correlation(&HMLFexh), correlation(&HMLFroh1), correlation(&HMLFroh100), correlation(&HMLFroh1000)); 
	fprintf(fsumR, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f ",
		correlation(&HMLmafHMLqtl), correlation(&HMLmafFhatI), correlation(&HMLmafFhatII), correlation(&HMLmafFhatIII), correlation(&HMLmafFhom), correlation(&HMLmafFexh), correlation(&HMLmafFroh1), correlation(&HMLmafFroh100), correlation(&HMLmafFroh1000)); 
	fprintf(fsumR, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f ",
		correlation(&HMLqtlFhatI), correlation(&HMLqtlFhatII), correlation(&HMLqtlFhatIII), correlation(&HMLqtlFhom), correlation(&HMLqtlFexh), correlation(&HMLqtlFroh1), correlation(&HMLqtlFroh100), correlation(&HMLqtlFroh1000)); 
	fprintf(fsumR, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f ",
		correlation(&FhatIFhatII), correlation(&FhatIFhatIII), correlation(&FhatIFhom), correlation(&FhatIFexh), correlation(&FhatIFroh1), correlation(&FhatIFroh100), correlation(&FhatIFroh1000)); 
	fprintf(fsumR, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f ",
		correlation(&FhatIIFhatIII), correlation(&FhatIIFhom), correlation(&FhatIIFexh), correlation(&FhatIIFroh1), correlation(&FhatIIFroh100), correlation(&FhatIIFroh1000)); 
	fprintf(fsumR, "%6.4f %6.4f %6.4f %6.4f %6.4f ",
		correlation(&FhatIIIFhom), correlation(&FhatIIIFexh), correlation(&FhatIIIFroh1), correlation(&FhatIIIFroh100), correlation(&FhatIIIFroh1000)); 
	fprintf(fsumR, "%6.4f %6.4f %6.4f %6.4f ",
		correlation(&FhomFexh), correlation(&FhomFroh1), correlation(&FhomFroh100), correlation(&FhomFroh1000)); 
	fprintf(fsumR, "%6.4f %6.4f %6.4f ",
		correlation(&FexhFroh1), correlation(&FexhFroh100), correlation(&FexhFroh1000)); 
	fprintf(fsumR, "%6.4f %6.4f ", correlation(&Froh1Froh100), correlation(&Froh1Froh1000)); 
	fprintf(fsumR, "%6.4f\n", correlation(&Froh100Froh1000));

	return(0);
}

/* ********************************************************************************************* */

lookfortext(s)
char *s;
{
   int len, i, curchar;
   char c;

   curchar = 0;
   len = 0;

   for (i=0; i<=100; i++)
   {
      if (s[i] == '\0') break;
      len++;
   }
   do
   {
      c = getc(fdat);

      if (c==s[curchar])
      {
         curchar++;
         if (curchar==len) return(0);
      }
      else curchar = 0;
   }
   while (c != EOF);
}

/* ********************************************************************************************* */

