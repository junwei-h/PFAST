/* This is program is to implement Transmission Travel time Tomography using Adjoint state method and FSM algorithm, 
*  based on the papers:
* [1] Zhao Hongkai, 2004, A Fast Sweeping Method for Eikonal Equaitons, Mathemtatics of Computiation, 
*     Vol. 74, No. 250, pp. 603-627  
* [2] Leung, Shingyu and Qian, Jianliang, 2006, An Adjoint State Method for Three-dimensional Transmission Traveltime
*	  Tomograhy using Frist-arrivals, Comm. Math. Sci., Vol. 4, No. 1, pp. 249-266.
*
* First created by Junwei Huang @ NRCan, Mar 08 2010, 
* PLEASE DO NOT DISTRIBUTE. PLEASE REFER OTHER PEOPLE TO THE AUTHOR. 

* Author	
* Dr. Jun-Wei Huang, 
* Geological Survey of Canada,
* 615 Booth Street, Ottawa, Ontario,
* Canada
* ph: +1 613 992 4968

* mailto:juhuang@nrcan.gc.ca
* Homepage: http://www.junweihuang.info

* If you want to publish results calculated with this program please
* give a reference to the following paper:
	Huang, J.W., and Bellefleur, G., 2011, Joint Transmission and Reflection Traveltime Tomography in Imaging Permafrost and Thermokarst lakes in Northwest Territories, Canada, The Leading Edge.
	Huang, J.W., and Bellefleur, G., 2012, Joint Transmission and Reflection Traveltime Tomography using the Fast Sweeping Method and the Adjoint-state Technique, Geophysical Journal International.
* ----------------------------------------------------------------

* History of important modifications:

* 21.02.2010	Version 1.0	original implementation of fast sweeping
				Jun-Wei Huang
* 08.03.2010	Implemented inversion  
				Jun-Wei Huang
* 22.03.2010    Implemented Parallel implementation
				Jun-Wei Huang
* 07.04.2010	Version 1.01 Implement Nonlinear Conjugate Gradient Method, Hestenes-Stiefel scheme
				Jun-Wei Huang
* 11.04.2010	Version 1.02 Implement Strong Wolfe Condition for line searching
				Jun-Wei Huang
* 12.04.2010	Version 1.03 Implement Nonlinear Conjugate Gradient Method, CGDECENT scheme according to
				"Harger, W.W. and Zhang, H.C., 2005, 
				A new Conjugate Gradient Method with Guaranteed Descent and an Efficient Line Search, 
				SIAM, J. OPTIM., Vo. 16, No. 1, pp. 170-192"
				Jun-Wei Huang
* 13.04.2010	Version 1.04 Implement L-BFGS Quasi-Newton Method
				Jun-Wei Huang
* 23.04.2010	version 1.05 Implement practical source and receiver geometry handling,see irr_fsm3d.h, irr_adjoint3d.h
				Jun-Wei Huang
*/
//#include "nr3.h"
#include "T_pfsm3d.h"
#include "T_padjoint3d.h"
#include <mpi.h>

/*--------FFT Accuracy Control Switch----------*/
#define FFT_NOFLOAT
#undef REAL
#ifdef TEST_FLOAT
# define REAL	float
# define fftn	fftnf
#else
# define REAL	double
#endif

using namespace std;

int main(int argc, char *argv[])
{	
	FILE *fp;
	int NX,NY,NZ,M=0,LM,schemetag=3,lsrc=3,NumofSrc,NumofRec,NumofObs,NumofR=0,count=0,Iterk=0,MAXIter=100,MAXLineSch=2,MAXSwp=100,ConvFlag=0,isresume=0;
	int *src,*srnum,*rec,sNX,sNY,rNX,rNY,ch0,ch1,chr;
	char *inpf,line[100],inp[20][100],tempar[8][80],vp0file[100],srcfile[100],recfile[100],obsfile[100],vpfile[100],lmdfile[100],p_vpfile[100],cmd[100];
	double dx,dy,dz,*T_obs,*nxyz,*T_cal,c1=1.0e-4,c2=0.9;
	double nux,nuy,nuz;
	float temp; 
	
	//===========start multi-processes======== 
	int my_rank,np;
	double Telapse,tic,lstime;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&np);
	
	MPI_Status status;
	if (my_rank==0)	Telapse=MPI_Wtime();
	if (my_rank==0)	info();
	if (argc==2){
		inpf=argv[1];LM=5;
		if (my_rank==0)	printf("Only Applicable to Optimization Method: L-BFGS\n\t...Default Memory Number m=%d\n",LM);
	}
	else if (argc==3){
		inpf=argv[1];LM=atoi(argv[2]);
		if (my_rank==0)	printf("Only Applicable to Optimization Method: L-BFGS\n\t...Specified Memory Number m=%d\n",LM);
	}
	else{
		if (my_rank==0)	 usage();
		return 1;
	}
//	printf("input arguments %d\ninfile: %s\noutfile: %s",argc,inpf,outf);

    if( (fp  = fopen( inpf, "r" )) == NULL ) 
		{if (my_rank==0) printf( "The file %s was not opened\n",inpf );}
   else
		{if (my_rank==0) printf( "Read input parameters from file %s\nUser Parameters are as follows:\n",inpf );}

	while(fgets(line,100,fp)){		
		//printf("%s",line);
		if(line[0]!='#'){
			strcpy(inp[M],line);
			M=M+1;			
		}
	}
	fclose(fp);

	///////////////////////////
	//Read Input parameters
	sscanf(inp[0],"%s",vp0file);
	
	sscanf(inp[1],"%s",srcfile);
	
	sscanf(inp[2],"%s%s",tempar,tempar+1);sNX=atoi(tempar[0]);sNY=atoi(tempar[1]);
	
	sscanf(inp[3],"%s",recfile);
	
	sscanf(inp[4],"%s%s",tempar,tempar+1);rNX=atoi(tempar[0]);rNY=atoi(tempar[1]);
	
	sscanf(inp[5],"%s%s%s%s%s%s",tempar,tempar+1,tempar+2,tempar+3,tempar+4,tempar+5);
	NX=atoi(tempar[0]);NY=atoi(tempar[1]);NZ=atoi(tempar[2]);dx=atof(tempar[3]);dy=atof(tempar[4]);dz=atof(tempar[5]);
	
	sscanf(inp[6],"%s",vpfile);
	
	sscanf(inp[7],"%s%s%s",tempar,tempar+1,tempar+2);
	nux=atof(tempar[0]);nuy=atof(tempar[1]);nuz=atof(tempar[2]);
	
	sscanf(inp[8],"%s%s",tempar,tempar+1);schemetag=atoi(tempar[0]);lsrc=atoi(tempar[1]); 
	//////////////////////////
	
	/////////////////////////
	//Report Input
	if (my_rank==0) printf("\nThe initial velocity file... %s\n",vp0file);
	if (vp0file[strlen(vp0file)-1]=='p') {
		isresume=0;
		if (my_rank==0) printf("|==========This is a NEW run===========|\n");
	}
	else {
		isresume=vp0file[strlen(vp0file)-1]-48;
		if (my_rank==0) printf("|==========This is a RESUME run===========|\n");
	}
	if (my_rank==0) printf("The source location file... %s\n",srcfile);
	if (my_rank==0) printf("The source binary file dimension... %d x %d\n",sNX,sNY);
	if (my_rank==0) printf("The receiver location file... %s\n",recfile);
	if (my_rank==0) printf("The receiver binary file dimension... %d x %d\n",rNX,rNY);
	if (my_rank==0) printf("The size of the model... %d X %d X %d\nSample Interval... %.2lfm x %.2lfm x %.2lfm\n",NX,NY,NZ,dx,dy,dz);
	if (my_rank==0)	printf("Final Velocity Model Output to %s\n",vpfile);
	if (my_rank==0)	printf("Regularization Parameters, (nux, nuy, nuz)=(%f, %f, %f)\n",nux,nuy,nuz);
	///////////////////////
	
#ifndef NoInverse //if undefined, this is an inverse run
	if (my_rank==0) {
		switch(schemetag){
			case 3:
				printf("L-BFGS quasi-newton method IN USE\n");
			break;
				
			case 1:
				printf("HS Nonlinear Conjugate Method IN USE\n");
			break;
				
			case 2:
				printf("CGDESCENT Nonlinear Conjugate Method IN USE\n");
			break;
			
			case 0:
				printf("Steepest Decent IN USE\n");
			break;
		}
		switch(lsrc){
			case 1:
				printf("Line Search Method: Secant method with exact Conditions\n");
			break;
				
			case 2:
				printf("Line Search Method: Secant method with Strong Wolfe Conditions\n");
			break;
				
			case 3:
				printf("Line Search Method: Cubic Interpolation method with Strong Wolfe Conditions\n");
			break;
			
			case 4:
				printf("Line Search Method: No Line Search, i.e. Steepest Decent\n");
			break;
		}
	}
#endif	
	double *vp=new (nothrow)double [NX*NY*NZ];//NX=number of rows, NY=number of colums
	double *T=new (nothrow)double [NX*NY*NZ];//store travel time
	double *lambda_i;//store single shot adjoint state
	double *lambda=new (nothrow)double [NX*NY*NZ];//store local shots adjoint state			
	double *srcArray=new (nothrow) double [sNX*sNY];
	double *recArray=new (nothrow) double [rNX*rNY];
	if (!T && !vp && !lambda) {printf("PE_%d: Not enough memory...abort...\n",my_rank);return 1;}
	////////////////////////////
	//To Load the velocity model 
	if (my_rank==0) {
		readmod(vp, vp0file, NX, NY,NZ);
		readmod(srcArray,srcfile,sNX,sNY,1);
		readmod(recArray,recfile,rNX,rNY,1);
	}
	MPI_Barrier( MPI_COMM_WORLD );
	MPI_Bcast( &vp[0],  NX*NY*NZ,  MPI_DOUBLE,  0,  MPI_COMM_WORLD);
	MPI_Bcast( &srcArray[0],  sNX*sNY,  MPI_DOUBLE,  0,  MPI_COMM_WORLD);
	MPI_Bcast( &recArray[0],  rNX*rNY,  MPI_DOUBLE,  0,  MPI_COMM_WORLD);
	///////////////////////////	
	
	///////////////////////////
	//Read NumofSrc
	NumofSrc=sNX;	
	src=new (nothrow)int [NumofSrc*DIM];
	srnum=new (nothrow)int [NumofSrc];
	
	int rmd=NumofSrc % np; //remainder
	int loc_nofs=(NumofSrc-(rmd))/np; //shots number per PE
	int *loc_srcid= new int [np];
	if (rmd!=0) { 
		if (my_rank==0) printf("!!!WARNING...Number of Shots %d is not divisible to %d PEs\n\n",NumofSrc,np);
		for (int i=0;i<np;i++){	
			if (i<rmd)
				loc_srcid[i]=(i+1)*(loc_nofs+1);
			else
				loc_srcid[i]=loc_srcid[i-1]+loc_nofs;					
		}
	}
	else {
		for (int i=0;i<np;i++)
			loc_srcid[i]=(i+1)*(loc_nofs);
	}
	
	printf("PE_%d (out of %d) handles %d shots: %d ~ %d\n",my_rank,np,loc_srcid[my_rank]-(my_rank-1>=0?loc_srcid[my_rank-1]:0),(my_rank-1>=0?loc_srcid[my_rank-1]:0)+1,loc_srcid[my_rank]);
	MPI_Barrier( MPI_COMM_WORLD );
	///////////////////////////
	//Read Source Locations
	if (my_rank==0) {
	
		printf("...reading source file:\n\t%s\n",srcfile);
		printf(" ----------------------------------\n");
		printf("| Total Number of Sources: %d\n",NumofSrc);
		printf("| index\trow\tcolumn\tnz\tNumofChn\n"); 
	
		for (int ii=0;ii<NumofSrc;ii++){
			src[ii*DIM+0]=(int)ceil(srcArray[ii*sNY+0]/dy)-1;//C++ start index from 0
			src[ii*DIM+1]=(int)ceil(srcArray[ii*sNY+1]/dx)-1;
			src[ii*DIM+2]=(int)ceil(srcArray[ii*sNY+2]/dz)-1;
			if (vp[src[ii*DIM+2]*NX*NY+src[ii*DIM+0]*NY+src[ii*DIM+1]]<400) {
			    printf("Source %d:\t%d\t%d\t%d is not in the model volume, please check...\n",ii,src[ii*DIM]+1,src[ii*DIM+1]+1,src[ii*DIM+2]+1);
				exit(1);
			}
			srnum[ii]=(int)ceil(srcArray[ii*sNY+3]);//channel number for this shot
			if (ii<=20 )
				printf("| %d\t%d\t%d\t%d\t%d\n",ii,src[ii*DIM]+1,src[ii*DIM+1]+1,src[ii*DIM+2]+1,srnum[ii]);
			else if (ii==21)
				printf("| ......\n");
			else if (ii==(NumofSrc-1) )
				printf("| %d\t%d\t%d\t%d\t%d\n",ii,src[ii*DIM]+1,src[ii*DIM+1]+1,src[ii*DIM+2]+1,srnum[ii]);
		}
		printf("| The End of Source Locations\n");
		printf(" ---------------------------------\n");
	}
	MPI_Barrier( MPI_COMM_WORLD );
	MPI_Bcast( &src[0], NumofSrc*DIM,  MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast( &srnum[0], NumofSrc,  MPI_INT, 0, MPI_COMM_WORLD);
	
	ch1=0;for (int j=(my_rank-1>=0?loc_srcid[my_rank-1]:0);j<loc_srcid[my_rank];j++) ch1+=srnum[j];
	//printf("PE_%d: ch1=%d\n",my_rank,ch1);
	T_cal=new (nothrow) double [ch1*(NumofR+1)];
	//T_cal=new (nothrow)double [NumofRec];
	if (!T_cal) {printf("PE_%d: Not enough memory for T_cal...sadly abort...\n",my_rank);return 1;}
	///////////////////////////////////
	
	///////////////////////////
	//Read Receivers
	NumofRec=rNX;
	rec=new (nothrow)int [NumofRec*DIM]; nxyz=new (nothrow)double [NumofRec*DIM];T_obs=new (nothrow)double [NumofRec];
	
	if (my_rank==0) {
		printf("...reading receiver file:\n\t%s\n",recfile);
		printf(" ----------------------------------\n");
		printf("| Total Number of Channels: %d\n",NumofRec);
		printf("| index\trow\tcolumn\tnz\t\tnormal_ny\tnormal_nx\tnormal_nz\tFirst Arrival\n"); 		
		
		for (int ii=0;ii<NumofRec;ii++){
			rec[ii*DIM]=(int)ceil(recArray[ii*rNY+0]/dy)-1;//C++ start index from 0
			rec[ii*DIM+1]=(int)ceil(recArray[ii*rNY+1]/dx)-1;
			rec[ii*DIM+2]=(int)ceil(recArray[ii*rNY+2]/dz)-1;			
			if (vp[rec[ii*DIM+2]*NX*NY+rec[ii*DIM+0]*NY+rec[ii*DIM+1]]<400) {
			    printf("Receiver %d:\t%d\t%d\t%d is not in the model volume, please check...\n",ii,rec[ii*DIM]+1,rec[ii*DIM+1]+1,rec[ii*DIM+2]+1);
				exit(1);
			}
			nxyz[ii*DIM]=recArray[ii*rNY+3];//C++ start index from 0
			nxyz[ii*DIM+1]=recArray[ii*rNY+4];
			nxyz[ii*DIM+2]=recArray[ii*rNY+5];
			
			T_obs[ii]=recArray[ii*rNY+6];
			
			if (ii<=20)
				printf("| %d\t%d\t%d\t%d\t\t%6.3f\t%6.3f\t%6.3f\t%6.3f\n",ii,rec[ii*DIM]+1,rec[ii*DIM+1]+1,rec[ii*DIM+2]+1,nxyz[ii*DIM],nxyz[ii*DIM+1],nxyz[ii*DIM+2],T_obs[ii]);
			else if (ii==21)
				printf("| ......\n");
			else if (ii==(NumofRec-1))
				printf("| %d\t%d\t%d\t%d\t\t%6.3f\t%6.3f\t%6.3f\t%6.3f\n",ii,rec[ii*DIM]+1,rec[ii*DIM+1]+1,rec[ii*DIM+2]+1,nxyz[ii*DIM],nxyz[ii*DIM+1],nxyz[ii*DIM+2],T_obs[ii]);
		}
		printf("| The End of Receiver Locations\n");
		printf(" ---------------------------------\n");
	}
	MPI_Barrier( MPI_COMM_WORLD );
	MPI_Bcast( &rec[0], NumofRec*DIM,  MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast( &nxyz[0], NumofRec*DIM,  MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//MPI_Bcast( &rectag[0], NX*NY*NZ,  MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast( &T_obs[0], NumofRec,  MPI_DOUBLE, 0, MPI_COMM_WORLD);
	delete [] recArray;recArray=NULL;
	delete [] srcArray;srcArray=NULL;
	
	strcpy(p_vpfile,"");
	if (isresume==0)
		strncat(p_vpfile,vp0file,strlen(vp0file)-3);//printf("%s\n",p_vpfile);
	else if (isresume<=9)
		strncat(p_vpfile,vp0file,strlen(vp0file)-4);
	else if (isresume<=99)
		strncat(p_vpfile,vp0file,strlen(vp0file)-5);
		
	#ifndef NoInverse //if undefined, this is an inverse run
	
	//////////////////////////
	//Start Iterating
	bool convflag=false;
	double alpha0=0.0,alpha1=5.0e3,alphak,f0,f1,beta;
	double *Djk1=new (nothrow)double [NX*NY*NZ],*Djk0=new (nothrow)double [NX*NY*NZ], *dk=new (nothrow) double [NX*NY*NZ],*Dj1=new (nothrow)double [NX*NY*NZ],*loc_vp=new (nothrow)double [NX*NY*NZ];
	double *sk=new (nothrow) double [NX*NY*NZ*LM], *yk=new (nothrow) double [NX*NY*NZ*LM],*rhoi=new (nothrow) double [LM],*q=new (nothrow) double [NX*NY*NZ],*ai=new (nothrow) double [LM],tmp,gamma;
	
	if (sk && yk){for (int i=0;i<NX*NY*NZ*LM;i++) {sk[i]=0.0;yk[i]=0.0;}}
	else {printf("PE%d: Memory Shortage for sk and yk...abort...\n",my_rank);return 1;}
	
	//For regularization	
	double global_E[MAXIter],loc_E=0.0;
	FILE *fe;
	
	if (isresume) {
		if (my_rank==0){ //for RESUME run				
			sprintf(lmdfile,"%s_bf_regu_%d.vp",p_vpfile,isresume);readmod(Djk1,lmdfile,NX,NY,NZ);
		}
		MPI_Bcast( &Djk1[0],  NX*NY*NZ,  MPI_DOUBLE,  0,  MPI_COMM_WORLD);
		for (int ii=0;ii<NX*NY*NZ;ii++)
				dk[ii]=Djk1[ii];
		Iterk=Iterk+isresume;
		if (my_rank==0) printf("Iterating No. %d resumed: (nux,nuy,nuz)=(%f,%f,%f)\n",Iterk,nux,nuy,nuz);
	}
//	else
//		strcpy(p_vpfile,vp0file);
	
	while (!convflag  && Iterk++<MAXIter){
		//nux=2.*nux/double(Iterk+1);nuy=2.*nuy/double(Iterk+1);nuz=2.*nuz/double(Iterk+1);
		if (my_rank==0) printf("Iterating No. %d started: (nux,nuy,nuz)=(%f,%f,%f)\n",Iterk,nux,nuy,nuz);
		tic=MPI_Wtime();		
		for (int kk=0;kk<NZ*NY*NX;kk++)	{lambda[kk]=0.0;Djk0[kk]=Djk1[kk];}
		loc_E=0.0;
		
		if (isresume) {int i=0;
			if (my_rank==0){
				sprintf(lmdfile,"%s_bf_regu_%d.vp",p_vpfile,Iterk);readmod(Djk1,lmdfile,NX,NY,NZ);
				if ((fe=fopen("global_E.tmp","r"))==NULL) printf(" The file global_E.tmp cannot be opened\n");
				i=0;while(fgets(line,100,fe)){
					sscanf(line,"%s",tempar);global_E[i++]=atof(tempar[0]);//printf("%s, %g\n",line,global_E[i]);i++;
				}
				fclose(fe);
				if ((fe=fopen("alpha1.tmp","r"))==NULL) printf(" The file alpha1.tmp cannot be opened\n");
				while(fgets(line,100,fe)){
					sscanf(line,"%s",tempar);alpha1=atof(tempar[0]);//printf("%g\n",alpha1);
				}
				fclose(fe);
				printf("Average Time Residual: %g ms\n",1000*sqrt(2*global_E[Iterk-1]/NumofRec));
				//printf("i=%d\n",i);
			}
			MPI_Bcast( &Djk1[0],  NX*NY*NZ,  MPI_DOUBLE,  0,  MPI_COMM_WORLD);
			MPI_Bcast( &global_E[0], i,  MPI_DOUBLE, 0, MPI_COMM_WORLD);			
			isresume=0;
		}
		else {
			for (int srcid=(my_rank-1>=0?loc_srcid[my_rank-1]:0);srcid<loc_srcid[my_rank];srcid++){
				fsm3d(vp,T,src,srcid,MAXSwp,dx,dy,dz,NX,NY,NZ,my_rank);
				lambda_i=new (nothrow)double [NX*NY*NZ];
				if (!lambda_i) {printf("PE_%d: Not enough memory for lambda_i...abort...\n",my_rank);return 1;}
				loc_E+=adjoint_fsm3d(lambda_i,T,T_obs,rec,nxyz,src,srnum,srcid,MAXSwp,dx,dy,dz,NX,NY,NZ,my_rank);

				for (int kk=0;kk<NX*NY*NZ;kk++) lambda[kk]=lambda[kk]-lambda_i[kk]/(vp[kk]*vp[kk]*vp[kk]);
							//lambda[kk]=lambda[kk]+lambda_i[kk];
				delete [] lambda_i;lambda_i=NULL;			
			}		
			MPI_Barrier( MPI_COMM_WORLD );
			MPI_Allreduce( &lambda[0], &Djk1[0], NX*NY*NZ, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
			MPI_Allreduce( &loc_E, &global_E[Iterk-1], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
			
			if (my_rank==0) {
				printf("Average Time Residual: %g ms\n",1000*sqrt(2*global_E[Iterk-1]/NumofRec));
				//fprintf(fe,"%f\n",global_E[Iterk-1]);
			}
			
			if (my_rank==0 ) {sprintf(lmdfile,"%s_bf_regu_%d.vp",p_vpfile,Iterk);writemod(Djk1,lmdfile,NX,NY,NZ);}
			//printf("PE%d: before regu3d\n",my_rank);
			
			regu3d(Djk1,vp,nux, nuy, nuz,dx,dy,dz, NX,NY,NZ);
			//printf("PE%d: after regu3d\n",my_rank);
			
			if (schemetag==3) {
				if (Iterk==1) alpha1=5e2*magnitudeOf(Djk1,NX,NY,NZ);
			}
			else alpha1=5e2*magnitudeOf(Djk1,NX,NY,NZ);
			//printf("PE%d: before magnitudeOf\n",my_rank);
			if (my_rank==0 ) {sprintf(lmdfile,"%s_af_regu_%d.vp",p_vpfile,Iterk); writemod(Djk1,lmdfile,NX,NY,NZ);}
			if (my_rank==0) {sprintf(lmdfile,"%s_GlobalErr.txt",p_vpfile); if ((fe=fopen(lmdfile,"w"))==NULL) printf(" The file %s cannot be opened\n",lmdfile);
				for (int i=0;i<Iterk;i++)
					fprintf(fe,"%g\n",global_E[i]);
				fclose(fe);}
		}
		/*
		1-HS Nonlinear Conjugate Method
		2-CGDESCENT Nonlinear Conjugate Method
		3-L-BFGS quasi-newton method
		*/
		schemetag=3;
		if (Iterk==1) {
			for (int ii=0;ii<NX*NY*NZ;ii++)
				dk[ii]=Djk1[ii];			
		}
		else {
			switch(schemetag)
			{
				case 1:
					beta=beta_HS(Djk1,Djk0,dk,NX,NY,NZ);
					for (int ii=0;ii<NX*NY*NZ;ii++)
						dk[ii]=Djk1[ii]+dk[ii]*beta;
					c1=1.0e-4;c2=0.1;
				break;
					
				case 2:
					beta=beta_CGDESCENT(Djk1,Djk0,dk,NX,NY,NZ);
					for (int ii=0;ii<NX*NY*NZ;ii++)
						dk[ii]=Djk1[ii]+dk[ii]*beta;
					c1=1.0e-4;c2=0.1;
				break;
				
				case 3:					
					//gamma=gamma_LBFGS(Djk1,Djk0,dk,alpha1,NX,NY,NZ);
					//q=new (nothrow) double [NX*NY*NZ];
					//if (q) {for (int j=0;j<NX*NY*NZ;j++) q[j]=Djk1[j];}
					//else {printf("PE_%d: Memory shortage for q...abort...\n",my_rank);return 1;}
					//ai=new (nothrow) double [Iterk-1-((Iterk-LM-1)<=0?0:Iterk-LM-1)];
					//if (!ai) {printf("PE_%d: Memory shortage for q...abort...\n",my_rank);return 1;}
					//vector <double> q(0.0,NX*NY*NZ);
					//vector <double> ai(0.0,Iterk-1-(Iterk-LM-1<=0?0:Iterk-LM-1));
					for (int j=0;j<NX*NY*NZ;j++) q[j]=Djk1[j];
					//printf("PE_%d: size of ai=%d\n",my_rank,Iterk-1-((Iterk-LM-1)<=0?0:Iterk-LM-1));
					gamma=gd2sk(Djk1,Djk0,dk,alpha1,rhoi,sk,yk,Iterk,LM,NX,NY,NZ);
					//if (Iterk-2>=LM-1) printf("PE_%d: gamma=%f, after gd2sk\n",my_rank,gamma);				
					//L_BFGS two-loop recursion to update searching direction dk
					for (int i=Iterk-2;i>=((Iterk-LM-1>=0)?Iterk-LM-1:0);i--)
					{
						ai[i%LM]=0.0;for (int j=0;j<NX*NY*NZ;j++) ai[i%LM]+=rhoi[i%LM]*sk[(i%LM)*NX*NY*NZ+j]*q[j];
						for (int j=0;j<NX*NY*NZ;j++) q[j]-=ai[i]*yk[(i%LM)*NX*NY*NZ+j];
					}
					for (int j=0;j<NX*NY*NZ;j++) dk[j]=gamma*q[j];
					//if (Iterk-2>=LM-1) printf("PE_%d: break 1\n",my_rank);
					//delete [] q;q=NULL;
					for (int i=((Iterk-LM-1>=0)?Iterk-LM-1:0);i<=Iterk-2;i++)
					{
						tmp=0.0;for (int j=0;j<NX*NY*NZ;j++) tmp+=rhoi[i%LM]*yk[(i%LM)*NX*NY*NZ+j]*dk[j];
						for (int j=0;j<NX*NY*NZ;j++) dk[j]+=sk[(i%LM)*NX*NY*NZ+j]*(ai[i%LM]-tmp);
					}
					//for (int j=0;j<NX*NY*NZ;j++) dk[j]=-dk[j];
					//if (Iterk-2>=LM-1) printf("PE_%d: break 2\n",my_rank);
					//delete [] ai;ai=NULL;
					
					alpha0=0.0;alpha1=1.0;
					c1=1.0e-4;c2=0.9;				
					
				break;
				
				/*case 0:
					for (int ii=0;ii<NX*NY*NZ;ii++)
						dk[ii]=Djk1[ii];
					lsrc=0;
				break;*/
				
				default:
					if (my_rank==0) printf("Numerical Optimization Method Unknown..abort...\n");
					return 1;
				break;
			}
		}
		//int lsrc=3;
		/*Line search method:
		1-Secant method with exact Conditions
		2-Line Search with Strong Wolfe Conditions
		3-Cubic Interpolation method, 1 iteration only
		0-No Line Search
		*/
		switch (lsrc)
		{
		case 1: 
			{lstime=MPI_Wtime();
			if (my_rank==0) printf("Secant Method Line Search for 'exact' condition (alpha1=%g)...\n",alpha1);
			//if (my_rank==0) printf("(alpha_0 alpha_1)=(%f %f)\n",alpha0,alpha1);
			count=0;alpha0=0.;loc_vp=new (nothrow) double [NX*NY*NZ];
			for (int i=0;i<NX*NY*NZ;i++) f0=f0+Djk1[i]*dk[i];
				
			while (fabs(alpha1-alpha0)/fabs(alpha0)>1e-2 && count++<MAXLineSch){			
				for (int i=0;i<NX*NY*NZ;i++){loc_vp[i]=vp[i]+alpha1*dk[i];lambda[i]=0.0;}
				for (int srcid=(my_rank-1>=0?loc_srcid[my_rank-1]:0);srcid<loc_srcid[my_rank];srcid++){
					fsm3d(loc_vp,T,src,srcid,MAXSwp,dx,dy,dz,NX,NY,NZ,my_rank);
					lambda_i=new (nothrow)double [NX*NY*NZ];
					if (!lambda_i) {printf("PE_%d: Not enough memory for lambda_i...abort...\n",my_rank);return 1;}
						//adjoint_fsm3d(lambda_i,T,T_obs,rec,nxyz,src,srcid,NumofSrc,NumofRec,MAXSwp,dx,dy,dz,NX,NY,NZ,my_rank);		
					adjoint_fsm3d(lambda_i,T,T_obs,rec,nxyz,src,srnum,srcid,MAXSwp,dx,dy,dz,NX,NY,NZ,my_rank);
					for (int ii=0;ii<NX*NY*NZ;ii++)	lambda[ii]=lambda[ii]-lambda_i[ii]/(loc_vp[ii]*loc_vp[ii]*loc_vp[ii]);	
					delete [] lambda_i;lambda_i=NULL;		
							
				}
				MPI_Barrier( MPI_COMM_WORLD );
				MPI_Allreduce( &lambda[0], &Dj1[0], NX*NY*NZ, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
				regu3d(Dj1,vp,nux, nuy, nuz,dx,dy,dz, NX,NY,NZ);
				//if (my_rank==0) {sprintf(lmdfile,"Dj1_after_regu%d.tmp",Iterk); writemod(Dj1,lmdfile,NX,NY,NZ);}
				f0=0.0;f1=0.0;
				for (int i=0;i<NX*NY*NZ;i++) f1=f1+Dj1[i]*dk[i];
				
				alphak=alpha1-(alpha1-alpha0)*f1/(f1-f0);alpha0=alpha1;alpha1=alphak;f0=f1;
				//if (my_rank==0) printf("f1=%20f, f0=%20f\n",f1,f0);
				
				//if (my_rank==0) printf("Updated (alpha_0 alpha_1)=(%f %f), changed by %f\n",alpha0,alpha1,fabs((alpha1-alpha0)/alpha0));
			}
			if (count>MAXLineSch) {
				if(my_rank==0) printf("Line Search Stopped, takes %g seconds, Stepsize=%f\n ...Updating velocity now...\n",MPI_Wtime()-lstime,alpha1);}
			else {
				if(my_rank==0) printf("Line Search Completed After %d Iterations, takes %g seconds, Stepsize=%f\n...Updating velocity now...\n",count,MPI_Wtime()-lstime,alpha1);}
		
			delete [] loc_vp;loc_vp=NULL;			
			break;}
		
		case 2: //with strong wolfe conditions
			{lstime=MPI_Wtime();			
			if (my_rank==0) printf("Secant Method Line Search for Strong Wolfe Condition (alpha1=%g)...\n",alpha1);
			
			bool isWolfe=false;
			double phi0=global_E[Iterk-1],phi1;
			double dphi0=0.0,dphi1=0.0;
			count=0;alpha0=0.0;
			//if (my_rank==0) printf("(alpha_0 alpha_1)=(%f %f)\n",alpha0,alpha1);
			dphi0=0.0; for (int i=0;i<NX*NY*NZ;i++) dphi0+=-dk[i]*Djk1[i];
			while (!isWolfe && count++<MAXLineSch ){	
				for (int i=0;i<NX*NY*NZ;i++) {loc_vp[i]=vp[i]+alpha1*dk[i];lambda[i]=0.0;}
				loc_E=0.0;
				for (int srcid=(my_rank-1>=0?loc_srcid[my_rank-1]:0);srcid<loc_srcid[my_rank];srcid++){
					fsm3d(loc_vp,T,src,srcid,MAXSwp,dx,dy,dz,NX,NY,NZ,my_rank);
					lambda_i=new (nothrow)double [NX*NY*NZ];
					if (!lambda_i) {printf("PE_%d: Not enough memory for lambda_i...abort...\n",my_rank);return 1;}
						//loc_E+=adjoint_fsm3d(lambda_i,T,T_obs,rec,nxyz,rectag,src,srcid,NumofSrc,NumofRec,MAXSwp,dx,dy,dz,NX,NY,NZ,my_rank);		
					loc_E+=adjoint_fsm3d(lambda_i,T,T_obs,rec,nxyz,src,srnum,srcid,MAXSwp,dx,dy,dz,NX,NY,NZ,my_rank);
					for (int ii=0;ii<NX*NY*NZ;ii++)lambda[ii]=lambda[ii]-lambda_i[ii]/(loc_vp[ii]*loc_vp[ii]*loc_vp[ii]);
					delete [] lambda_i;lambda_i=NULL;
					
				}
				MPI_Barrier( MPI_COMM_WORLD );
				MPI_Allreduce( &lambda[0], &Dj1[0], NX*NY*NZ, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
				MPI_Allreduce( &loc_E, &phi1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
				regu3d(Dj1,vp,nux, nuy, nuz,dx,dy,dz, NX,NY,NZ);
				dphi1=0.0;for (int i=0;i<NX*NY*NZ;i++) dphi1+=-dk[i]*Dj1[i];
				//if (count==1) {alpha1=alpha1-(alpha1-alpha0)*dphi1/(dphi1-dphi0);}
				if(my_rank==0) printf("(alpha0,alpha1,phi0,phi1,dphi0,dphi1)=(%g %g %g %g %g %g\n",alpha0,alpha1,phi0,phi1,dphi0,dphi1);
				
				if (phi1<=phi0+c1*alpha1*dphi0 && fabs(dphi1)<=-c2*dphi0)
					isWolfe=true;
				else
					alpha1=alpha1-(alpha1-alpha0)*dphi1/(dphi1-dphi0);
				//if (my_rank==0) printf("(alpha_0 alpha_1)=(%g %g)\n",alpha0,alpha1);
			}
			if (my_rank==0) {if ((fe=fopen("alpha1.tmp","w"))==NULL) printf(" The file alpha1.tmp cannot be opened\n");
				fprintf(fe,"%g\n",alpha1);fclose(fe);}
			if (count>MAXLineSch) {
				if(my_rank==0) printf("Line Search Stopped, takes %g seconds, Stepsize=%f\n ...Updating velocity now...\n",MPI_Wtime()-lstime,alpha1);
			}
			else {
				if(my_rank==0) printf("Wolfe condition satisfied After %d Iterations, takes %g seconds, Stepsize=%f\n...Updating velocity now...\n",count,MPI_Wtime()-lstime,alpha1);
				//if(my_rank==0) printf("(alpha0,alpha1,phi0,phi1,dphi0,dphi1)=(%g %g %g %g %g %g\n",alpha0,alpha1,phi0,phi1,dphi0,dphi1);
			}
		
			//delete [] loc_vp;loc_vp=NULL;
			break;}
		
		case 3: //cubic interpolation
		{
			lstime=MPI_Wtime();			
			if (my_rank==0) printf("Cubic Interpolation Line Searching (alpha1=%g)...\n",alpha1);

			bool isWolfe=false;
			double phi0=global_E[Iterk-1],phi1;
			double dphi0=0.0,dphi1=0.0;
			count=0;alpha0=0.0;
			dphi0=0.0; for (int i=0;i<NX*NY*NZ;i++) dphi0+=-dk[i]*Djk1[i];
			while (!isWolfe && count++<MAXLineSch ){
				for (int i=0;i<NX*NY*NZ;i++) {loc_vp[i]=vp[i]+alpha1*dk[i];lambda[i]=0.0;}			
				loc_E=0.0;
				for (int srcid=(my_rank-1>=0?loc_srcid[my_rank-1]:0);srcid<loc_srcid[my_rank];srcid++){
					fsm3d(loc_vp,T,src,srcid,MAXSwp,dx,dy,dz,NX,NY,NZ,my_rank);
					lambda_i=new (nothrow)double [NX*NY*NZ];
					if (!lambda_i) {printf("PE_%d: Not enough memory for lambda_i...abort...\n",my_rank);return 1;}
					//loc_E+=adjoint_fsm3d(lambda_i,T,T_obs,rec,nxyz,rectag,src,srcid,NumofSrc,NumofRec,MAXSwp,dx,dy,dz,NX,NY,NZ,my_rank);
					loc_E+=adjoint_fsm3d(lambda_i,T,T_obs,rec,nxyz,src,srnum,srcid,MAXSwp,dx,dy,dz,NX,NY,NZ,my_rank);		
					for (int ii=0;ii<NX*NY*NZ;ii++)	lambda[ii]=lambda[ii]-lambda_i[ii]/(loc_vp[ii]*loc_vp[ii]*loc_vp[ii]);
					delete [] lambda_i;lambda_i=NULL;
				}
				MPI_Barrier( MPI_COMM_WORLD );
				MPI_Allreduce( &lambda[0], &Dj1[0], NX*NY*NZ, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
				MPI_Allreduce( &loc_E, &phi1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
				regu3d(Dj1,vp,nux, nuy, nuz,dx,dy,dz, NX,NY,NZ);
				dphi1=0.0;for (int i=0;i<NX*NY*NZ;i++) dphi1+=-dk[i]*Dj1[i];				
				if (phi1<=phi0+c1*alpha1*dphi0 && fabs(dphi1)<=-c2*dphi0)
					isWolfe=true;
				else
					alpha1=cubicInterp(alpha0,alpha1,phi0,phi1,dphi0,dphi1);

				if(my_rank==0) printf("(alpha0,alpha1,phi0,phi1,dphi0,dphi1)=(%g %g %g %g %g %g\n",alpha0,alpha1,phi0,phi1,dphi0,dphi1);
				//alpha1=zoom(alpha0,alpha1,phi0,phi1,dphi0,dphi1);
			}
			if (count>MAXLineSch) {
				if(my_rank==0) printf("Line Search Stopped, takes %g seconds, Stepsize=%f\n ...Updating velocity now...\n",MPI_Wtime()-lstime,alpha1);
			}
			else {
				if(my_rank==0) printf("Wolfe condition satisfied After %d Iterations, takes %g seconds, Stepsize=%f\n...Updating velocity now...\n",count,MPI_Wtime()-lstime,alpha1);
			}
			break;}
		case 0:
			if (my_rank==0) printf("No Line Searching, Stepsize=%f\n",alpha1);
			break;
		default:
			if (my_rank==0) printf("Line search tag unknown..abort...\n");
			return 1;
			break;
		}
		
		updatevp(vp,dk,alpha1,NX,NY,NZ);		
		if (my_rank==0) { 
			sprintf(lmdfile,"%s.vp%d",p_vpfile,Iterk);
			writemod(vp,lmdfile,NX,NY,NZ);
		}
		convflag=check_conv(dk,alpha1,vp,1e-5,global_E,Iterk,NX,NY,NZ);
		MPI_Barrier( MPI_COMM_WORLD );
		if (my_rank==0) printf("--------Interating No. %d takes %10.4f seconds---------------\n\n",Iterk,MPI_Wtime()-tic);
	}	
	
	MPI_Barrier( MPI_COMM_WORLD );
	if(my_rank==0) {
		writemod(vp, vpfile, NX, NY,NZ);
		sprintf(lmdfile,"%s_AveRes_T.txt",p_vpfile);
		if ((fe=fopen(lmdfile,"w"))==NULL)
			printf(" The file %s cannot be opened\n",lmdfile);
		for (int i=0;i<Iterk;i++)
			fprintf(fe,"%g\n",1000*sqrt(2*global_E[i]/NumofRec));
		fclose(fe);
	
	printf("Total Elapse Time: %10.4f seconds (%8.4f hours) \n",MPI_Wtime()-Telapse, (MPI_Wtime()-Telapse)/3600.0);}

#else
// Forward modeling and save the traveltime at receiver locations

	for (int srcid=(my_rank-1>=0?loc_srcid[my_rank-1]:0);srcid<loc_srcid[my_rank];srcid++){
		ch0=0;for (int j=(my_rank-1>=0?loc_srcid[my_rank-1]:0);j<srcid;j++) ch0+=srnum[j];
		chr=0;for (int j=0;j<srcid;j++) chr+=srnum[j];
		//printf("PE_%d:ch0=%d,chr=%d\n",my_rank,ch0,chr);
		fsm3d(vp,T,src,srcid,MAXSwp,dx,dy,dz,NX,NY,NZ,my_rank);		
		for (int idx=0;idx<srnum[srcid];idx++) {	
			int i=rec[(idx+chr)*DIM+0],j=rec[(idx+chr)*DIM+1],k=rec[(idx+chr)*DIM+2];
			//int i=rec[(idx+chr)*DIM+0],j=rec[(idx+chr)*DIM+1];
			//if (TT_obs[(ch0+idx)*(NumofR+1)+0]<100)
				//T_cal[(ch0+idx)*(NumofR+1)+0]=T[i*NY+j];
				T_cal[(ch0+idx)*(NumofR+1)+0]=T[k*NX*NY+i*NY+j];
			//else
				//T_cal[(ch0+idx)*(NumofR+1)+0]=INF;
		}
	}
	sprintf(lmdfile,"Tr3D.%d",my_rank);writemod(T_cal, lmdfile, ch1, NumofR+1,1);
	MPI_Barrier( MPI_COMM_WORLD );
	
	if(my_rank==0) {
		for (int i=0;i<np;i++){
			if (i==0) {sprintf(cmd,"cat Tr3D.%d > %s_Tr3D.bin ",i,p_vpfile);system(cmd);}
			else {sprintf(cmd,"cat Tr3D.%d >> %s_Tr3D.bin ",i,p_vpfile);system(cmd);}
			sprintf(cmd,"rm Tr3D.%d",i);system(cmd);		
		}
		printf("Total Elapse Time: %10.4f seconds (%8.4f hours) \n",MPI_Wtime()-Telapse, (MPI_Wtime()-Telapse)/3600.0);
	}

#endif
	MPI_Finalize();
	return 0;

}