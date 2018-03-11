/* This is program has no output file but computation time!
*  The purpose of this program is to test the scalibity with multiple processors. 
*  For application to real cases, please use other programs in this package. 
*  
*  --A message from the Author of this program
*  First created by Junwei Huang @ NRCan, March 28 2011
*/

#include "RT_pfsm2d.h"
#include "RT_padjoint2d.h"
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
	int NX,NY,M=0,LM,schemetag=1,RTtag=3,lsrc,NumofSrc,NumofRec,NumofR,count=0,Iterk=0,MAXIter=500,MAXLineSch=5,MAXSwp=500,ConvFlag=0;
	int *src,*srnum,*rec,sNX,sNY,rNX,rNY,gmNX,gmNY,ch0,ch1,chr;
	char *inpf,line[100],inp[20][100],tempar[8][80],vp0file[100],srcfile[100],recfile[100],refgeom[100],obsfile[100],vpfile[100],twfile[100],lmdfile[100],p_obsfile[100],cmd[100];
	double dx,dz,*nxy,*TT_obs,*TR_obs, c1=1.0e-4,c2=0.9;//*WR;
	double nux=80.,nuy=50.;
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
#ifndef NoInverse	
	if (argc==2){
		inpf=argv[1];LM=5;RTtag=3;
		if (my_rank==0)	printf("Optimization Method: L-BFGS\nDefault Memory Number m=%d (Ignore if L-BFGS is not used)\n",LM);
	}
	else if (argc==3){
		inpf=argv[1];LM=atoi(argv[2]);RTtag=3;
		if (my_rank==0)	printf("Optimization Method: L-BFGS\nSpecified Memory Number m=%d (Ignore if L-BFGS is not used)\n",LM);
	}
	else if (argc==4){
		inpf=argv[1];LM=atoi(argv[2]);RTtag=atoi(argv[3]);
		if (my_rank==0)	printf("Optimization Method: L-BFGS\nSpecified Memory Number m=%d (Ignore if L-BFGS is not used)\n",LM);
	}
	else{
		if (my_rank==0)	 usage();
		return 1;
	}
	
	if (my_rank==0) {
		printf(" *****************************************************************\n");
		switch(RTtag){
			case 3:
				printf("Joint Inversion of Reflection and Transmission\n");
			break;
				
			case 1:
				printf("Inversion of Transmission ONLY\n");
			break;
				
			case 2:
				printf("Inversion of Reflection ONLY\n");
			break;
		}
		printf(" *****************************************************************\n");
	}
#else
	inpf=argv[1];
#endif
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
	
	sscanf(inp[5],"%s",refgeom);
	
	sscanf(inp[6],"%s%s",tempar,tempar+1);gmNX=atoi(tempar[0]);gmNY=atoi(tempar[1]);
	
	sscanf(inp[7],"%s%s%s%s",tempar,tempar+1,tempar+2,tempar+3);
	NX=atoi(tempar[0]);NY=atoi(tempar[1]);dx=atof(tempar[2]);dz=atof(tempar[3]);
	
	sscanf(inp[8],"%s",vpfile);
	
	sscanf(inp[9],"%s%s",tempar,tempar+1);	
	nux=atof(tempar[0]);nuy=atof(tempar[1]);
	
	sscanf(inp[10],"%s",twfile);
	
	sscanf(inp[11],"%s%s",tempar,tempar+1);schemetag=atoi(tempar[0]);lsrc=atoi(tempar[1]); 
	//////////////////////////
		
	///////////////////////
	//Report Input
	if (my_rank==0) 
	{
		printf("\nThe initial velocity file... %s\n",vp0file);
		printf("The source location file... %s\n",srcfile);
		printf("The receiver location file... %s\n",recfile);
		printf("The source binary file dimension... %d x %d\n",sNX,sNY);
		printf("The receiver binary file dimension... %d x %d\n",rNX,rNY);
		printf("The reflector geometry file... %s\n",refgeom);
		printf("The reflector geometry binary file dimension... %d x %d\n",gmNX,gmNY);
		printf("The size of the model... %d X %d\nSample Interval... %.2lfm x %.2lfm\n",NX,NY,dx,dz);
		printf("Final Velocity Model Output to %s\n",vpfile);
		printf("The weighting factor matrix file... %s\n\n",twfile);
	}
	//////////////////////
	
	if (my_rank==0)	{
		switch(schemetag){
			case 1:
				printf("HS Nonlinear Conjugate Method In Use\n");
			break;
				
			case 2:
				printf("CGDESCENT Nonlinear Conjugate Method In Use\n");
			break;
				
			case 3:
				printf("L-BFGS quasi-newton method In Use\n");
			break;
			
			case 0: 
				printf("Steepest Decent Method In Use\n");
			break;
		}
		switch(lsrc){
			case 1:
				printf("Secant method with exact Conditions\n");
			break;
				
			case 2:
				printf("Line Search with Strong Wolfe Conditions\n");
			break;
				
			case 3:
				printf("Cubic Interpolation method with Strong Wolfe Conditions\n");
			break;
			
			case 0: 
				printf("No Line Search for steepest decent only\n");
			break;
		}
	}
	
	
	NumofR=(int)ceil((double)gmNX/NY);					//number of reflectors
	double *vp=new (nothrow)double [NX*NY];				//NX=number of rows, NY=number of colums
	double *TT=new (nothrow)double [NX*NY];				//store Transmission Travel time in the whole model
	double *TTd=new (nothrow)double [NX*NY*NumofR];		//store downgoing travel time
	double *TTu=new (nothrow)double [NX*NY*NumofR];		//store upgoing travel time
	double *Tr;											//store travel time at receiver locations, not identical among processors
	double *Z=new (nothrow)double [gmNX*gmNY]; 			//store reflector geometry
	double *lambdaT_i,*lambdaR_i;						//store single shot adjoint state of transmission and reflection wave
	double *lambdaT,*lambdaR,*lambdaRT;					//store local shots adjoint state within one processing unit of transmission and reflection wave
	double *srcArray=new (nothrow) double [sNX*sNY];	//store source z, x, number of channel
	double *recArray=new (nothrow) double [rNX*rNY];	//store receiver z, x, -nz, -nx, FB, R1, R2,...Rn
		 
	if (!TT && !TTd && !TTu && !vp && !Z && !srcArray && !recArray ) {
		printf("PE_%d: Not enough memory...sadly abort...\n",my_rank);
		return 1;
	}
	
	////////////////////////////
	//To Load the velocity model 
	if (my_rank==0) {
		readmod(vp, vp0file, NX, NY);
		readmod(srcArray,srcfile,sNX,sNY);
		readmod(recArray,recfile,rNX,rNY);
		readmod(Z,refgeom,gmNX,gmNY);
	}
	MPI_Barrier( MPI_COMM_WORLD );
	MPI_Bcast( &vp[0],  NX*NY,  MPI_DOUBLE,  0,  MPI_COMM_WORLD);
	MPI_Bcast( &srcArray[0],  sNX*sNY,  MPI_DOUBLE,  0,  MPI_COMM_WORLD);
	MPI_Bcast( &recArray[0],  rNX*rNY,  MPI_DOUBLE,  0,  MPI_COMM_WORLD);
	MPI_Bcast( &Z[0],  gmNX*gmNY,  MPI_DOUBLE,  0,  MPI_COMM_WORLD);
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
		printf("| index\trow\tcolumn\tNumofChn\n"); 
	
		for (int ii=0;ii<NumofSrc;ii++){
			src[ii*DIM+0]=(int)ceil(srcArray[ii*sNY+0]/dz)-1;//C++ start index from 0
			src[ii*DIM+1]=(int)ceil(srcArray[ii*sNY+1]/dx)-1;
			srnum[ii]=(int)(srcArray[ii*sNY+2]);
			if (ii<=20 )
				printf("| %d\t%d\t%d\t%d\n",ii,src[ii*DIM]+1,src[ii*DIM+1]+1,srnum[ii]);
			else if (ii==21)
				printf("| ......\n");
			else if (ii==(NumofSrc-1) )
				printf("| %d\t%d\t%d\t%d\n",ii,src[ii*DIM]+1,src[ii*DIM+1]+1,srnum[ii]);
		}
		printf("| The End of Source Locations\n");
		printf(" ---------------------------------\n");
	}
	MPI_Barrier( MPI_COMM_WORLD );
	MPI_Bcast( &src[0], NumofSrc*DIM,  MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast( &srnum[0], NumofSrc,  MPI_INT, 0, MPI_COMM_WORLD);
	
	ch1=0;for (int j=(my_rank-1>=0?loc_srcid[my_rank-1]:0);j<loc_srcid[my_rank];j++) ch1+=srnum[j];
	//printf("PE_%d: ch1=%d\n",my_rank,ch1);
	Tr=new (nothrow) double [ch1*(NumofR+1)];
	if (!Tr) {printf("PE_%d: Not enough memory for Tr...sadly abort...\n",my_rank);return 1;}
	///////////////////////////////////
	
	///////////////////////////
	//Read Receivers
	NumofRec=rNX;
	rec=new (nothrow)int [NumofRec*DIM]; nxy=new (nothrow)double [NumofRec*DIM];
	TT_obs=new (nothrow)double [NumofRec];TR_obs=new (nothrow)double [NumofRec*NumofR];
	
	if (my_rank==0) {
		printf("...reading receiver file:\n\t%s\n",recfile);
		printf(" ----------------------------------\n");
		printf("| Total Number of Channels: %d\n",NumofRec);
		printf("| index\trow\tcolumn\t\tnormal_ny\tnormal_nx\tTT\tTR1......TRn\n"); 		
		
		for (int ii=0;ii<NumofRec;ii++){
			rec[ii*DIM]=(int)ceil(recArray[ii*rNY+0]/dz)-1;//C++ start index from 0
			rec[ii*DIM+1]=(int)ceil(recArray[ii*rNY+1]/dx)-1;		
				
			nxy[ii*DIM]=recArray[ii*rNY+2];//C++ start index from 0
			nxy[ii*DIM+1]=recArray[ii*rNY+3];
			
			TT_obs[ii]=recArray[ii*rNY+4];
			for (int it=0;it<NumofR;it++) TR_obs[ii*NumofR+it]=recArray[ii*rNY+5+it];
			
			if (ii<=20)
				printf("| %d\t%d\t%d\t\t%6.3f\t%6.3f\t%6.3f\t%6.3f......%6.3f\n",ii,rec[ii*DIM]+1,rec[ii*DIM+1]+1,nxy[ii*DIM],nxy[ii*DIM+1],TT_obs[ii],TR_obs[ii*NumofR+0],TR_obs[ii*NumofR+NumofR-1]);
			else if (ii==21)
				printf("| ......\n");
			else if (ii==(NumofRec-1))
				printf("| %d\t%d\t%d\t\t%6.3f\t%6.3f\t%6.3f\t%6.3f......%6.3f\n",ii,rec[ii*DIM]+1,rec[ii*DIM+1]+1,nxy[ii*DIM],nxy[ii*DIM+1],TT_obs[ii],TR_obs[ii*NumofR+0],TR_obs[ii*NumofR+NumofR-1]);
		}
		printf("| The End of Receiver Locations\n");
		printf(" ---------------------------------\n");
	}
	MPI_Barrier( MPI_COMM_WORLD );
	MPI_Bcast( &rec[0], NumofRec*DIM,  MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast( &nxy[0], NumofRec*DIM,  MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//MPI_Bcast( &rectag[0], NX*NY,  MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast( &TT_obs[0], NumofRec,  MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast( &TR_obs[0], NumofRec*NumofR,  MPI_DOUBLE, 0, MPI_COMM_WORLD);
	delete [] recArray;recArray=NULL;
	delete [] srcArray;srcArray=NULL;
	
	///////////////////////////
	//Read Reflector Geometry	
	if (my_rank==0) {
		printf("...reading reflector geometry file:\n\t%s\n",refgeom);
		printf(" ----------------------------------\n");
		printf("| Total Number of reflector geometry points: %d\n",gmNX);
		printf("| Total Number of reflector layers: %d\n",NumofR);
		printf("| index\tdz\t\tref_ny\tref_nx\n"); 		
		
		for (int ii=0;ii<gmNX;ii++){
		
			if (ii<=20)
				printf("| %d\t%6.3f\t%6.3f\t%6.3f\n",ii,Z[ii*3+0],Z[ii*3+1],Z[ii*3+2]);
			else if (ii==21)
				printf("| ......\n");
			else if (ii==(gmNX-1))
				printf("| %d\t%6.3f\t%6.3f\t%6.3f\n",ii,Z[ii*3+0],Z[ii*3+1],Z[ii*3+2]);
		}
		printf("| The End of Reflector Geometries\n");
		printf(" ---------------------------------\n");
	}
	MPI_Barrier( MPI_COMM_WORLD );
	
	#ifndef NoInverse //if undefined, this is an inverse run
	
	//////////////////////////
	//Start Iterating
	bool convflag=false;
	double alpha0=0.0,alpha1=5.0e3,alphak,f0,f1,beta;
	double *Djk1=new (nothrow)double [NX*NY],*Djk0=new (nothrow)double [NX*NY], *dk=new (nothrow) double [NX*NY],*Dj1=new (nothrow)double [NX*NY],*loc_vp=new (nothrow)double [NX*NY];
	double *sk=new (nothrow) double [NX*NY*LM], *yk=new (nothrow) double [NX*NY*LM],*rhoi=new (nothrow) double [LM],*q=new (nothrow) double [NX*NY],*ai=new (nothrow) double [LM],tmp,gamma;
	double *WR=new (nothrow)double [NX*NY];	//store spatial varying weighting factors, NX=number of rows, NY=number of colums
	
	////////////////////////////
	//To Load the weighting factor model 
	if (my_rank==0) {
		readmod(WR, twfile, NX, NY);
	}
	MPI_Barrier( MPI_COMM_WORLD );
	MPI_Bcast( &WR[0],  NX*NY,  MPI_DOUBLE,  0,  MPI_COMM_WORLD);
	
	if (sk && yk){for (int i=0;i<NX*NY*LM;i++) {sk[i]=0.0;yk[i]=0.0;}}
	else {printf("PE%d: Memory Shortage for sk and yk...sadly abort...\n",my_rank);return 1;}
	
	//For regularization
	double global_E[MAXIter],global_Et[MAXIter],global_Er[MAXIter],loc_E=0.0,loc_Er=0.0,loc_Et=0.0;
	FILE *fe;
	
	while (!convflag  && Iterk++<MAXIter){
		//nux=2.*nux/double(Iterk+1);nuy=2.*nuy/double(Iterk+1);
		if (my_rank==0) printf("Iterating No. %d: (nux,nuy)=(%f,%f)\n",Iterk,nux,nuy);
		tic=MPI_Wtime();
		lambdaR=new (nothrow)double [NX*NY];lambdaT=new (nothrow)double [NX*NY];lambdaRT=new (nothrow)double [NX*NY];
		for (int kk=0;kk<NY*NX;kk++)	{lambdaR[kk]=0.0;lambdaT[kk]=0.0;lambdaRT[kk]=0.0;Djk0[kk]=Djk1[kk];}
		loc_E=0.0;loc_Er=0.0;loc_Et=0.0;
		
		for (int srcid=(my_rank-1>=0?loc_srcid[my_rank-1]:0);srcid<loc_srcid[my_rank];srcid++){
			//ch0=0;for (int j=(my_rank-1>=0?loc_srcid[my_rank-1]:0);j<srcid;j++) ch0+=srnum[j];
			RT_fsm2d(vp, TT, TTd, TTu, Z, src,rec, srnum, srcid, NumofR, MAXSwp, dx,dz, NX,NY,my_rank);
			lambdaT_i=new (nothrow)double [NX*NY];
			loc_Et+=Tadjoint_fsm2d(lambdaT_i,TT,TT_obs,rec,nxy,src,srnum,srcid,MAXSwp,dx,dz,NX,NY,my_rank);
			for (int kk=0;kk<NX*NY;kk++) lambdaT[kk]=lambdaT[kk]-lambdaT_i[kk]/(vp[kk]*vp[kk]*vp[kk]);
			delete [] lambdaT_i;lambdaT_i=NULL;
			
			lambdaR_i=new (nothrow)double [NX*NY];
			loc_Er+=Radjoint_fsm2d(lambdaR_i,TTd,TTu, TR_obs,rec,nxy,src,srnum,srcid,NumofR,MAXSwp,dx,dz,NX,NY,my_rank);
			for (int kk=0;kk<NX*NY;kk++) lambdaR[kk]=lambdaR[kk]-lambdaR_i[kk]/(vp[kk]*vp[kk]*vp[kk]);
			delete [] lambdaR_i;lambdaR_i=NULL;	

			//RTweight(WR,lambdaRT,lambdaR,lambdaT,&loc_E,&loc_Er,&loc_Et,NX,NY);
			for (int kk=0;kk<NX*NY;kk++) lambdaRT[kk]=lambdaT[kk]*(1.-WR[kk])+lambdaR[kk]*WR[kk];
			loc_E=loc_Er+loc_Et;
		}		
		MPI_Barrier( MPI_COMM_WORLD );
		if (RTtag==3) {//Joint
			MPI_Allreduce( &lambdaRT[0], &Djk1[0], NX*NY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );			
			}
		else if (RTtag==1) {//Transmission
			MPI_Allreduce( &lambdaT[0], &Djk1[0], NX*NY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );					
			}
		else if (RTtag==2) {//Reflection
			MPI_Allreduce( &lambdaR[0], &Djk1[0], NX*NY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );			
			}
		else if (my_rank==0)
			{printf("RTtag is unrecognized!...abort...\n");exit(1);}
			
		MPI_Allreduce( &loc_E, &global_E[Iterk-1], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
		MPI_Allreduce( &loc_Et, &global_Et[Iterk-1], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
		MPI_Allreduce( &loc_Er, &global_Er[Iterk-1], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
		
		double *gRtmp=new double [NX*NY];
		double *gTtmp=new double [NX*NY];
		MPI_Allreduce( &lambdaR[0], &gRtmp[0], NX*NY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
		MPI_Allreduce( &lambdaT[0], &gTtmp[0], NX*NY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

		if (my_rank==0) {
			printf("Average Time Residual (TR): %g ms\t",1000*sqrt(2*global_E[Iterk-1]/NumofRec));
			printf("(T): %g ms\t",1000*sqrt(2*global_Et[Iterk-1]/NumofRec));
			printf("(R): %g ms\n",1000*sqrt(2*global_Er[Iterk-1]/NumofRec));
			//fprintf(fe,"%f\n",global_E[Iterk-1]);
		}
		/*!!!!!!writemod commands are all commented out!!!!!!!!*/
		//if (my_rank==0 ) {sprintf(lmdfile,"%s_bf_regu_RT%d",vp0file,Iterk);writemod(Djk1,lmdfile,NX,NY);}
		regu2d(Djk1,vp,nux, nuy, dx,dz, NX,NY); 
		if (schemetag==3) {
			if (Iterk==1) alpha1=5e2*magnitudeOf(Djk1,NX,NY);
		}
		else alpha1=5e2*magnitudeOf(Djk1,NX,NY);
		
		//if (my_rank==0 ) {sprintf(lmdfile,"%s_af_regu_RT%d",vp0file,Iterk); writemod(Djk1,lmdfile,NX,NY);}
		
		//if (my_rank==0 ) {sprintf(lmdfile,"%s_bf_regu_R%d",vp0file,Iterk);writemod(gRtmp,lmdfile,NX,NY);}
		regu2d(gRtmp,vp,nux, nuy, dx,dz, NX,NY);
		//if (my_rank==0 ) {sprintf(lmdfile,"%s_af_regu_R%d",vp0file,Iterk); writemod(gRtmp,lmdfile,NX,NY);}
		
		//if (my_rank==0 ) {sprintf(lmdfile,"%s_bf_regu_T%d",vp0file,Iterk);writemod(gTtmp,lmdfile,NX,NY);}
		regu2d(gTtmp,vp,nux, nuy, dx,dz, NX,NY);
		//if (my_rank==0 ) {sprintf(lmdfile,"%s_af_regu_T%d",vp0file,Iterk); writemod(gTtmp,lmdfile,NX,NY);}
		delete [] gRtmp;gRtmp=NULL;	
		delete [] gTtmp;gTtmp=NULL;	
		if (my_rank==0) {
			sprintf(lmdfile,"%s_GlobalErr_TR_T_R.txt",vp0file); 
			if ((fe=fopen(lmdfile,"w"))==NULL) printf(" The file %s cannot be opened\n",lmdfile);
				for (int i=0;i<Iterk;i++)
					fprintf(fe,"%g\t%g\t%g\n",global_E[i],global_Et[i],global_Er[i]);
				fclose(fe);
		}
		//delete [] lambdaR;lambdaR=NULL;	
		//delete [] lambdaT;lambdaT=NULL;
		/*
		1-HS Nonlinear Conjugate Method
		2-CGDESCENT Nonlinear Conjugate Method
		3-L-BFGS quasi-newton method
		*/
		//schemetag=1;
		if (Iterk==1) {
			for (int ii=0;ii<NX*NY;ii++)
				dk[ii]=Djk1[ii];			
		}
		else {
			switch(schemetag)
			{
				case 1:
					beta=beta_HS(Djk1,Djk0,dk,NX,NY);
					for (int ii=0;ii<NX*NY;ii++)
						dk[ii]=Djk1[ii]+dk[ii]*beta;
					c1=1.0e-4;c2=0.1;
				break;
					
				case 2:
					beta=beta_CGDESCENT(Djk1,Djk0,dk,NX,NY);
					for (int ii=0;ii<NX*NY;ii++)
						dk[ii]=Djk1[ii]+dk[ii]*beta;
					c1=1.0e-4;c2=0.1;
				break;
				
				case 3:					
					for (int j=0;j<NX*NY;j++) q[j]=Djk1[j];
					//printf("PE_%d: size of ai=%d\n",my_rank,Iterk-1-((Iterk-LM-1)<=0?0:Iterk-LM-1));
					gamma=gd2sk(Djk1,Djk0,dk,alpha1,rhoi,sk,yk,Iterk,LM,NX,NY);
					//if (Iterk-2>=LM-1) printf("PE_%d: gamma=%f, after gd2sk\n",my_rank,gamma);
					//L_BFGS two-loop recursion to update searching direction dk
					for (int i=Iterk-2;i>=((Iterk-LM-1>=0)?Iterk-LM-1:0);i--)
					{
						ai[i%LM]=0.0;for (int j=0;j<NX*NY;j++) ai[i%LM]+=rhoi[i%LM]*sk[(i%LM)*NX*NY+j]*q[j];
						for (int j=0;j<NX*NY;j++) q[j]-=ai[i]*yk[(i%LM)*NX*NY+j];
					}
					for (int j=0;j<NX*NY;j++) dk[j]=gamma*q[j];
					//if (Iterk-2>=LM-1) printf("PE_%d: break 1\n",my_rank);
					//delete [] q;q=NULL;
					for (int i=((Iterk-LM-1>=0)?Iterk-LM-1:0);i<=Iterk-2;i++)
					{
						tmp=0.0;for (int j=0;j<NX*NY;j++) tmp+=rhoi[i%LM]*yk[(i%LM)*NX*NY+j]*dk[j];
						for (int j=0;j<NX*NY;j++) dk[j]+=sk[(i%LM)*NX*NY+j]*(ai[i%LM]-tmp);
					}
					//for (int j=0;j<NX*NY;j++) dk[j]=-dk[j];
					//if (Iterk-2>=LM-1) printf("PE_%d: break 2\n",my_rank);
					//delete [] ai;ai=NULL;
					
					alpha0=0.0;alpha1=1.0;
					c1=1.0e-4;c2=0.9;				
					
				break;
				
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
		3-Cubic Interpolation method with Strong Wolfe Conditions
		0-No Line Search
		*/
		switch (lsrc)
		{
		case 1: 
			{lstime=MPI_Wtime();
			if (my_rank==0) printf("Secant Method Line Search for 'exact' condition (alpha1=%g)...\n",alpha1);
			//if (my_rank==0) printf("(alpha_0 alpha_1)=(%f %f)\n",alpha0,alpha1);
			count=0;alpha0=0.e3;loc_vp=new (nothrow) double [NX*NY];
			for (int i=0;i<NX*NY;i++) f0=f0+Djk1[i]*dk[i];
				
			while (fabs(alpha1-alpha0)/fabs(alpha0)>1e-2 && count++<MAXLineSch){			
				for (int i=0;i<NX*NY;i++){loc_vp[i]=vp[i]+alpha1*dk[i];lambdaRT[i]=0.0;lambdaR[i]=0.0;lambdaT[i]=0.0;}
				loc_E=0.0;loc_Er=0.0;loc_Et=0.0;
				for (int srcid=(my_rank-1>=0?loc_srcid[my_rank-1]:0);srcid<loc_srcid[my_rank];srcid++){
					//ch0=0;for (int j=(my_rank-1>=0?loc_srcid[my_rank-1]:0);j<srcid;j++) ch0+=srnum[j];
					RT_fsm2d(loc_vp, TT, TTd, TTu, Z,src,rec, srnum,srcid, NumofR, MAXSwp, dx,dz, NX,NY,my_rank);
					lambdaT_i=new (nothrow)double [NX*NY];
					loc_Et+=Tadjoint_fsm2d(lambdaT_i,TT,TT_obs,rec,nxy,src,srnum,srcid,MAXSwp,dx,dz,NX,NY,my_rank);
					for (int kk=0;kk<NX*NY;kk++) lambdaT[kk]=lambdaT[kk]-lambdaT_i[kk]/(loc_vp[kk]*loc_vp[kk]*loc_vp[kk]);
					delete [] lambdaT_i;lambdaT_i=NULL;
					
					lambdaR_i=new (nothrow)double [NX*NY];
					loc_Er+=Radjoint_fsm2d(lambdaR_i,TTd,TTu,TR_obs,rec,nxy,src,srnum,srcid,NumofR,MAXSwp,dx,dz,NX,NY,my_rank);
					for (int kk=0;kk<NX*NY;kk++) lambdaR[kk]=lambdaR[kk]-lambdaR_i[kk]/(loc_vp[kk]*loc_vp[kk]*loc_vp[kk]);
					delete [] lambdaR_i;lambdaR_i=NULL;	
					
					//RTweight(WR,lambdaRT,lambdaR,lambdaT,&loc_E,&loc_Er,&loc_Et,NX,NY);
					for (int kk=0;kk<NX*NY;kk++) lambdaRT[kk]=lambdaT[kk]*(1.-WR[kk])+lambdaR[kk]*WR[kk];
					loc_E=loc_Er+loc_Et;
				}		
				MPI_Barrier( MPI_COMM_WORLD );
				if (RTtag==3) //Joint
					MPI_Allreduce( &lambdaRT[0], &Dj1[0], NX*NY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
				else if (RTtag==1) //Transmission
					MPI_Allreduce( &lambdaT[0], &Dj1[0], NX*NY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );			
				else if (RTtag==2) //Reflection
					MPI_Allreduce( &lambdaR[0], &Dj1[0], NX*NY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
				else if (my_rank==0)
					{printf("RTtag is unrecognized!...abort...\n");exit(1);}
					
				regu2d(Dj1,vp,nux, nuy, dx,dz, NX,NY);
				//if (my_rank==0) {sprintf(lmdfile,"Dj1_after_regu%d.tmp",Iterk); writemod(Dj1,lmdfile,NX,NY);}
				f0=0.0;f1=0.0;
				for (int i=0;i<NX*NY;i++) f1=f1+Dj1[i]*dk[i];
				
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
			if (my_rank==0) printf("Secant Method Line Search with Strong Wolfe Condition (alpha1=%g)...\n",alpha1);
			
			bool isWolfe=false;
			double phi0=global_E[Iterk-1],phi1;
			double dphi0=0.0,dphi1=0.0;
			count=0;alpha0=0.0;
			//if (my_rank==0) printf("(alpha_0 alpha_1)=(%f %f)\n",alpha0,alpha1);
			dphi0=0.0; for (int i=0;i<NX*NY;i++) dphi0+=-dk[i]*Djk1[i];
			while (!isWolfe && count++<MAXLineSch ){	
				for (int i=0;i<NX*NY;i++) {loc_vp[i]=vp[i]+alpha1*dk[i];lambdaRT[i]=0.0;lambdaR[i]=0.0;lambdaT[i]=0.0;}
				loc_E=0.0;loc_Er=0.0;loc_Et=0.0;
				for (int srcid=(my_rank-1>=0?loc_srcid[my_rank-1]:0);srcid<loc_srcid[my_rank];srcid++){
					//ch0=0;for (int j=(my_rank-1>=0?loc_srcid[my_rank-1]:0);j<srcid;j++) ch0+=srnum[j];
					RT_fsm2d(loc_vp, TT, TTd, TTu, Z,src,rec, srnum,srcid, NumofR, MAXSwp, dx,dz, NX,NY,my_rank);
					lambdaT_i=new (nothrow)double [NX*NY];
					loc_Et+=Tadjoint_fsm2d(lambdaT_i,TT,TT_obs,rec,nxy,src,srnum,srcid,MAXSwp,dx,dz,NX,NY,my_rank);
					for (int kk=0;kk<NX*NY;kk++) lambdaT[kk]=lambdaT[kk]-lambdaT_i[kk]/(loc_vp[kk]*loc_vp[kk]*loc_vp[kk]);
					delete [] lambdaT_i;lambdaT_i=NULL;
					
					lambdaR_i=new (nothrow)double [NX*NY];
					loc_Er+=Radjoint_fsm2d(lambdaR_i,TTd,TTu,TR_obs,rec,nxy,src,srnum,srcid,NumofR,MAXSwp,dx,dz,NX,NY,my_rank);
					for (int kk=0;kk<NX*NY;kk++) lambdaR[kk]=lambdaR[kk]-lambdaR_i[kk]/(loc_vp[kk]*loc_vp[kk]*loc_vp[kk]);
					delete [] lambdaR_i;lambdaR_i=NULL;	
					
					//RTweight(WR,lambdaRT,lambdaR,lambdaT,&loc_E,&loc_Er,&loc_Et,NX,NY);
					for (int kk=0;kk<NX*NY;kk++) lambdaRT[kk]=lambdaT[kk]*(1.-WR[kk])+lambdaR[kk]*WR[kk];
					loc_E=loc_Er+loc_Et;
				}
				MPI_Barrier( MPI_COMM_WORLD );
				if (RTtag==3) {//Joint
					MPI_Allreduce( &lambdaRT[0], &Dj1[0], NX*NY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
					MPI_Allreduce( &loc_E, &phi1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
				}
				else if (RTtag==1) {//Transmission
					MPI_Allreduce( &lambdaT[0], &Dj1[0], NX*NY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );		
					MPI_Allreduce( &loc_Et, &phi1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
				}
				else if (RTtag==2) {//Reflection
					MPI_Allreduce( &lambdaR[0], &Dj1[0], NX*NY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
					MPI_Allreduce( &loc_Er, &phi1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
				}
				else if (my_rank==0)
					{printf("RTtag is unrecognized!...abort...\n");exit(1);}
				//MPI_Allreduce( &lambdaRT[0], &Dj1[0], NX*NY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
				//MPI_Allreduce( &loc_E, &phi1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
				regu2d(Dj1,vp,nux, nuy, dx,dz, NX,NY);
				dphi1=0.0;for (int i=0;i<NX*NY;i++) dphi1+=-dk[i]*Dj1[i];
				//if (count==1) {alpha1=alpha1-(alpha1-alpha0)*dphi1/(dphi1-dphi0);}
				if(my_rank==0) printf("(alpha0,alpha1,phi0,phi1,dphi0,dphi1)=(%g %g %g %g %g %g\n",alpha0,alpha1,phi0,phi1,dphi0,dphi1);
				
				if (phi1<=phi0+c1*alpha1*dphi0 && fabs(dphi1)<=-c2*dphi0)
					isWolfe=true;
				else
					alpha1=alpha1-(alpha1-alpha0)*dphi1/(dphi1-dphi0);
				//if (my_rank==0) printf("(alpha_0 alpha_1)=(%g %g)\n",alpha0,alpha1);
			}
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
			if (my_rank==0) {
				if (schemetag==3) printf("Cubic Interpolation Line Search with Strong Wolfe Condition (alpha1=%g)...\n",alpha1);
				else printf("Cubic Interpolation Line Search (alpha1=%g)...\n",alpha1);
			}
			
			bool isWolfe=false;
			double phi0=global_E[Iterk-1],phi1;
			double dphi0=0.0,dphi1=0.0;
			count=0;alpha0=0.0;
			dphi0=0.0; for (int i=0;i<NX*NY;i++) dphi0+=-dk[i]*Djk1[i];
			if (schemetag==3){ //LBFGS, Strong Wolfe Condition Applied
				while (!isWolfe && count++<MAXLineSch ){
					for (int i=0;i<NX*NY;i++) {loc_vp[i]=vp[i]+alpha1*dk[i];lambdaRT[i]=0.0;lambdaR[i]=0.0;lambdaT[i]=0.0;}		
					loc_E=0.0;loc_Er=0.0;loc_Et=0.0;
					for (int srcid=(my_rank-1>=0?loc_srcid[my_rank-1]:0);srcid<loc_srcid[my_rank];srcid++){
							//ch0=0;for (int j=(my_rank-1>=0?loc_srcid[my_rank-1]:0);j<srcid;j++) ch0+=srnum[j];
						RT_fsm2d(loc_vp, TT, TTd, TTu, Z,src,rec, srnum,srcid, NumofR, MAXSwp, dx,dz, NX,NY,my_rank);
						lambdaT_i=new (nothrow)double [NX*NY];
						loc_Et+=Tadjoint_fsm2d(lambdaT_i,TT,TT_obs,rec,nxy,src,srnum,srcid,MAXSwp,dx,dz,NX,NY,my_rank);
						for (int kk=0;kk<NX*NY;kk++) lambdaT[kk]=lambdaT[kk]-lambdaT_i[kk]/(loc_vp[kk]*loc_vp[kk]*loc_vp[kk]);
						delete [] lambdaT_i;lambdaT_i=NULL;
						
						lambdaR_i=new (nothrow)double [NX*NY];
						loc_Er+=Radjoint_fsm2d(lambdaR_i,TTd,TTu,TR_obs,rec,nxy,src,srnum,srcid,NumofR,MAXSwp,dx,dz,NX,NY,my_rank);
						for (int kk=0;kk<NX*NY;kk++) lambdaR[kk]=lambdaR[kk]-lambdaR_i[kk]/(loc_vp[kk]*loc_vp[kk]*loc_vp[kk]);
						delete [] lambdaR_i;lambdaR_i=NULL;	
						
						//RTweight(WR,lambdaRT,lambdaR,lambdaT,&loc_E,&loc_Er,&loc_Et,NX,NY);
						for (int kk=0;kk<NX*NY;kk++) lambdaRT[kk]=lambdaT[kk]*(1.-WR[kk])+lambdaR[kk]*WR[kk];
						loc_E=loc_Er+loc_Et;
					}
					MPI_Barrier( MPI_COMM_WORLD );
					if (RTtag==3) {//Joint
							MPI_Allreduce( &lambdaRT[0], &Dj1[0], NX*NY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
							MPI_Allreduce( &loc_E, &phi1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
					}
						else if (RTtag==1) {//Transmission
							MPI_Allreduce( &lambdaT[0], &Dj1[0], NX*NY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );		
							MPI_Allreduce( &loc_Et, &phi1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
					}
						else if (RTtag==2) {//Reflection
							MPI_Allreduce( &lambdaR[0], &Dj1[0], NX*NY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
							MPI_Allreduce( &loc_Er, &phi1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
					}
					else if (my_rank==0)
						{printf("RTtag is unrecognized!...abort...\n");exit(1);}
		
					regu2d(Dj1,vp,nux, nuy, dx,dz, NX,NY);
					dphi1=0.0;for (int i=0;i<NX*NY;i++) dphi1+=-dk[i]*Dj1[i];
					
					if (phi1<=phi0+c1*alpha1*dphi0 && fabs(dphi1)<=-c2*dphi0)
						isWolfe=true;
					else
						alpha1=cubicInterp(alpha0,alpha1,phi0,phi1,dphi0,dphi1);

					if(my_rank==0) printf("(alpha0,alpha1,phi0,phi1,dphi0,dphi1)=(%g %g %g %g %g %g\n",alpha0,alpha1,phi0,phi1,dphi0,dphi1);
				}
				if (count>MAXLineSch) {
					if(my_rank==0) printf("Line Search Stopped, takes %g seconds, Stepsize=%f\n ...Updating velocity now...\n",MPI_Wtime()-lstime,alpha1);
				}
				else {
					if(my_rank==0) printf("Wolfe condition satisfied After %d Iterations, takes %g seconds, Stepsize=%f\n...Updating velocity now...\n",count,MPI_Wtime()-lstime,alpha1);
				}
				}
			else {
				for (int i=0;i<NX*NY;i++) {loc_vp[i]=vp[i]+alpha1*dk[i];lambdaRT[i]=0.0;lambdaR[i]=0.0;lambdaT[i]=0.0;}		
					loc_E=0.0;loc_Er=0.0;loc_Et=0.0;
					for (int srcid=(my_rank-1>=0?loc_srcid[my_rank-1]:0);srcid<loc_srcid[my_rank];srcid++){
							//ch0=0;for (int j=(my_rank-1>=0?loc_srcid[my_rank-1]:0);j<srcid;j++) ch0+=srnum[j];
						RT_fsm2d(loc_vp, TT, TTd, TTu, Z,src,rec, srnum,srcid, NumofR, MAXSwp, dx,dz, NX,NY,my_rank);
						lambdaT_i=new (nothrow)double [NX*NY];
						loc_Et+=Tadjoint_fsm2d(lambdaT_i,TT,TT_obs,rec,nxy,src,srnum,srcid,MAXSwp,dx,dz,NX,NY,my_rank);
						for (int kk=0;kk<NX*NY;kk++) lambdaT[kk]=lambdaT[kk]-lambdaT_i[kk]/(loc_vp[kk]*loc_vp[kk]*loc_vp[kk]);
						delete [] lambdaT_i;lambdaT_i=NULL;
						
						lambdaR_i=new (nothrow)double [NX*NY];
						loc_Er+=Radjoint_fsm2d(lambdaR_i,TTd,TTu,TR_obs,rec,nxy,src,srnum,srcid,NumofR,MAXSwp,dx,dz,NX,NY,my_rank);
						for (int kk=0;kk<NX*NY;kk++) lambdaR[kk]=lambdaR[kk]-lambdaR_i[kk]/(loc_vp[kk]*loc_vp[kk]*loc_vp[kk]);
						delete [] lambdaR_i;lambdaR_i=NULL;	
						
						//RTweight(WR,lambdaRT,lambdaR,lambdaT,&loc_E,&loc_Er,&loc_Et,NX,NY);
						for (int kk=0;kk<NX*NY;kk++) lambdaRT[kk]=lambdaT[kk]*(1.-WR[kk])+lambdaR[kk]*WR[kk];
						loc_E=loc_Er+loc_Et;
					}
					MPI_Barrier( MPI_COMM_WORLD );
					if (RTtag==3) {//Joint
							MPI_Allreduce( &lambdaRT[0], &Dj1[0], NX*NY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
							MPI_Allreduce( &loc_E, &phi1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
					}
						else if (RTtag==1) {//Transmission
							MPI_Allreduce( &lambdaT[0], &Dj1[0], NX*NY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );		
							MPI_Allreduce( &loc_Et, &phi1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
					}
						else if (RTtag==2) {//Reflection
							MPI_Allreduce( &lambdaR[0], &Dj1[0], NX*NY, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
							MPI_Allreduce( &loc_Er, &phi1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
					}
					else if (my_rank==0)
						{printf("RTtag is unrecognized!...abort...\n");exit(1);}
		
					regu2d(Dj1,vp,nux, nuy, dx,dz, NX,NY);
					dphi1=0.0;for (int i=0;i<NX*NY;i++) dphi1+=-dk[i]*Dj1[i];
					alpha1=cubicInterp(alpha0,alpha1,phi0,phi1,dphi0,dphi1);
					if(my_rank==0) printf("(alpha0,alpha1,phi0,phi1,dphi0,dphi1)=(%g %g %g %g %g %g\n",alpha0,alpha1,phi0,phi1,dphi0,dphi1);
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
		updatevp(vp,dk,alpha1,NX,NY);
		
		if (my_rank==0)  {
			sprintf(lmdfile,"%s_%d",vp0file,Iterk);
			//writemod(vp,lmdfile,NX,NY);
			//sprintf(lmdfile,"%s_%d.gTimeRes",vp0file,Iterk);
			//writemod(g_TimeRes, lmdfile, NumofRec, NumofR+1+DIM);
		}
		
		convflag=check_conv(dk,alpha1,vp,global_E,Iterk,10.0,NX,NY);
		MPI_Barrier( MPI_COMM_WORLD );
		if (my_rank==0) printf("--------Interating No. %d takes %10.4f seconds---------------\n\n",Iterk,MPI_Wtime()-tic);
	}	
	
	MPI_Barrier( MPI_COMM_WORLD );
	if(my_rank==0) {
		//writemod(vp, vpfile, NX, NY);
		sprintf(lmdfile,"%s_AveRes_TR_T_R.txt",vp0file);
		if ((fe=fopen(lmdfile,"w"))==NULL)
			printf(" The file Err.txt cannot be opened\n");
		for (int i=0;i<Iterk;i++)
			fprintf(fe,"%g\t%g\t%g\n",1000*sqrt(2*global_E[i]/NumofRec),1000*sqrt(2*global_Et[i]/NumofRec),1000*sqrt(2*global_Er[i]/NumofRec));
		fclose(fe);
		
		printf("Total Elapse Time: %10.4f seconds (%8.4f hours) \n",MPI_Wtime()-Telapse, (MPI_Wtime()-Telapse)/3600.0);
		RTinfo(RTtag);
	}
	#else
	for (int srcid=(my_rank-1>=0?loc_srcid[my_rank-1]:0);srcid<loc_srcid[my_rank];srcid++){
		ch0=0;for (int j=(my_rank-1>=0?loc_srcid[my_rank-1]:0);j<srcid;j++) ch0+=srnum[j];
		chr=0;for (int j=0;j<srcid;j++) chr+=srnum[j];
		//printf("PE_%d:ch0=%d,chr=%d\n",my_rank,ch0,chr);
		RT_fsm2d(vp, TT, TTd, TTu, Z,src,rec, srnum,srcid, NumofR, MAXSwp, dx,dz, NX,NY,my_rank);
		for (int idx=0;idx<srnum[srcid];idx++) {			
			int i=rec[(idx+chr)*DIM+0],j=rec[(idx+chr)*DIM+1];
			//if (TT_obs[(ch0+idx)*(NumofR+1)+0]<100)
				Tr[(ch0+idx)*(NumofR+1)+0]=TT[i*NY+j];
			//else
				//Tr[(ch0+idx)*(NumofR+1)+0]=INF;
		}
		for (int ir=0;ir<NumofR;ir++){ 
			for (int idx=0;idx<srnum[srcid];idx++) {			
				int i=rec[(idx+chr)*DIM+0],j=rec[(idx+chr)*DIM+1];
				//if (TR_obs[(ch0+idx)*(NumofR+1)+ir+1]>0)
					Tr[(ch0+idx)*(NumofR+1)+ir+1]=TTu[ir*NX*NY+i*NY+j];
				//else
				//	Tr[(ch0+idx)*(NumofR+1)+ir+1]=0;
			}
		}
		#ifndef NoShotGather
		//sprintf(lmdfile,"Src%d_TTd.2d",srcid);writemod(TTd,lmdfile,NumofR,NX*NY);
		//sprintf(lmdfile,"Src%d_TTu.2d",srcid);writemod(TTu,lmdfile,NumofR,NX*NY);
		#endif
	}
	//sprintf(lmdfile,"Tr.%d",my_rank);writemod(Tr, lmdfile, ch1, NumofR+1);
	MPI_Barrier( MPI_COMM_WORLD );
	
	if(my_rank==0) {
		/*for (int i=0;i<np;i++){
			if (i==0) {sprintf(cmd,"cat Tr.%d > Tr.bin ",i);system(cmd);}
			else {sprintf(cmd,"cat Tr.%d >> Tr.bin ",i);system(cmd);}
			sprintf(cmd,"rm Tr.%d",i);system(cmd);		
		}*/
		printf("Total Elapse Time: %10.4f seconds (%8.4f hours) \n",MPI_Wtime()-Telapse, (MPI_Wtime()-Telapse)/3600.0);}
	#endif
	return 0;
}
