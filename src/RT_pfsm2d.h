#ifndef _FIM2D_FUNCTIONS
#define _FIM2D_FUNCTIONS

#ifndef MAX
    #define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
    #define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif

#ifndef INF
    #define INF 9999999
#endif

#ifndef DIM
    #define DIM 2
#endif
/*
#ifndef ROUND
#define ROUND(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))
#endif
*/
#ifndef PI
	#define PI (3.141592653589793)
#endif

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <string.h>
//#include <algorithm>
//#include <vector>

using namespace std;

void readmod(double *mod, char *file, int nrows, int ncols)
{
	FILE *fp;
	if ((fp=fopen(file,"rb"))==NULL)
		printf(" The file %s cannot be opened\n",file);
	printf(" \n...reading file:\n\t %s\n",file);
	fread(mod, sizeof(double), nrows*ncols, fp);
	fclose(fp);
}

void writemod(double *mod, char *file, int nrows, int ncols)
{
	FILE *fp;
	if ((fp=fopen(file,"wb"))==NULL)
		printf(" The file %s cannot be opened\n",file);
	printf(" \n...writing file:\n\t %s\n",file);
	fwrite(mod, sizeof(double), nrows*ncols, fp);
	fclose(fp);
}

bool notin(int *node,int *CurrentNodes,int Num_CurrentNodes){
//to test if one node is in CurrentNodes

	clock_t tic=clock();
	for (int i=0;i<Num_CurrentNodes;i++)
		if (CurrentNodes[i*DIM+0]==node[0] && CurrentNodes[i*DIM+1]==node[1] && CurrentNodes[i*DIM+2]==node[2])
			{if ((double)(clock()-tic)/CLOCKS_PER_SEC > 1.0 ) printf("notin function takes %8.5f seconds\n", (double)(clock()-tic)/CLOCKS_PER_SEC);
			return false;}
	if ((double)(clock()-tic)/CLOCKS_PER_SEC > 1.0 )  printf("notin function takes %8.5f seconds\n", (double)(clock()-tic)/CLOCKS_PER_SEC);
	return true;		
}

double eikslv2d(double *U,double *vp,int *node,double dx,double dz,int nrows, int ncols)
{//to solve the time field U at node
	double a,b,q,ta,tb,tc,da=dx,db=dz,s=1.0/vp[node[0]*ncols+node[1]];

	if (node[0]-1<0)              //node is at the upper edge
		b=U[(node[0]+1)*ncols+node[1]];
	else if (node[0]+1>=nrows)    //node is at the lower edge
		b=U[(node[0]-1)*ncols+node[1]];
	else							//node is in the model
		b=MIN(U[(node[0]-1)*ncols+node[1]],U[(node[0]+1)*ncols+node[1]]);
	
	if (node[1]-1<0)              //node is at the left edge
		a=U[node[0]*ncols+node[1]+1];//printf("a=%f, take U(%d,%d+1)=%f\n",a,node[0],node[1],U[(node[0])*ncols+node[1]+1]);}
	else if (node[1]+1>=ncols)   //node is at the right edge
		a=U[node[0]*ncols+node[1]-1];//printf("a=%f, take U(%d,%d-1)=%f\n",a,node[0],node[1],U[(node[0])*ncols+node[1]-1]);}
	else                        //node is in the model
		a=MIN(U[node[0]*ncols+node[1]+1],U[node[0]*ncols+node[1]-1]);//printf("a=%f, take min(U(%d,%d+1)=%f,U(%d,%d-1)=%f\n",a, node[0],node[1],U[(node[0])*ncols+node[1]+1],node[0],node[1],U[(node[0])*ncols+node[1]-1]);}
	//if (isnan(a) || isnan(b)) {printf("(b,a)=(%f,%f) at node U(%d,%d)=%f\n",b,a,node[0],node[1],U[(node[0])*ncols+node[1]]);exit(1);}
	if (b>a) {swap(a,b);swap(da,db);} //to make sure a>>
	if (!(b<=a)) {printf("Error in sorting (b,a)=(%f,%f) at node U(%d,%d)=%f\n",b,a,node[0],node[1],U[(node[0])*ncols+node[1]]);exit(1);}
	
	//////////////////
	q=b+db*s;
	if (q<=a) return q;
	ta=1.0/da/da+1.0/db/db;tb=-2.*(a/da/da+b/db/db);tc=a*a/da/da+b*b/db/db-s*s;
	
	if ((tb*tb-4.*ta*tc)>=0) 
		q=(-tb+sqrt(tb*tb-4*ta*tc))/2./ta;
	else 
		q=INF;
		
	return q;
	/*if (fabs(a-b)>=pow(vp[node[0]*ncols+node[1]],-1)*simds)
		q=MIN(a,b)+pow(vp[node[0]*ncols+node[1]],-1)*simds;
	else
		q=(a+b+sqrt(2*pow(vp[node[0]*ncols+node[1]],-1)*simds*pow(vp[node[0]*ncols+node[1]],-1)*simds-(a-b)*(a-b)))/2.0;		
	return q;*/
}

bool isconv(double *p,double *U, int nrows, int ncols){
	for (int k=0;k<nrows*ncols;k++)
		if (fabs(p[k]-U[k])>1e-9) return false;			
	return true;
}

void RT_fsm2d(double *vp, double *U, double *TTd, double *TTu, double *Z, int *src, int *rec, int *srnum,int srcid, int NumofR, int MAXSwp, double dx,double dz, int nrows, int ncols,int pe)
{
	/*
	vp: velocity (NX x NY )
	U: transmission travel time through the whole model (NX x NY)
	TTd: Time field of downgoing wave, each reflector has one
	TTu: Time field of upgoing wave, each reflector has one
	Z: reflectors (NumofR*ncols*nrows x 3)=[dz,-ny,-nx]
	src: source locations (sy,sx)
	rec: receiver locations (ry,rx)
	srnum: number of receivers per shot (Numofshots x 1)
	srcid: source ID (scalor)
	NumofR: Number of reflectors (scalor)
	MAXSwp: maximum sweep time, default 500
	dx,dz: grid interval (scalor) along x and z directions
	nrows, ncols, nz: NX,NY,NZ
	pe: processor ID (scalor)
	*/
	
	//printf("|\tPE_%d Called Eikonal Equation Solver for source #%d (%d,%d,%d)...\n",pe,srcid+1,src[srcid*DIM+0]+1,src[srcid*DIM+1]+1,src[srcid*DIM+2]+1);
	int count=0,OneNode[DIM],k2,*k2vec=new int [ncols],chr=0;
	double *p=new double [nrows*ncols],qx,qz,dnx,dnz,rz,tmpU;
	FILE *fp;
	char filename[100];
	//double *U=new double [nrows*ncols];
	if (!p) {printf("Memory Shotage\n");exit (1);}
	clock_t tic=clock(),tmp;//start time recording
	
	for (int i=0;i<ncols;i++) k2vec[i]=-1;
	//////////////////////
	//Start Sweepting for reflections
	for (int ir=0;ir<NumofR;ir++) {
		tmp=clock();
		//if (pe==0) printf("Reflections from Layer: %d\n",ir+1);
		//============Initialize Time Field for Point Source==========
		for (int ii=0;ii<nrows;ii++){
			for (int jj=0;jj<ncols;jj++){
				if (ii==src[srcid*DIM+0] && jj==src[srcid*DIM+1]) 
					U[ii*ncols+jj]=0.0;
				else U[ii*ncols+jj]=INF;
				p[ii*ncols+jj]=0.0;
			}
		}
	
		//////////////////////
		//Start Sweepting
		count=0;
		while (!isconv(p,U,nrows,ncols) && count++<MAXSwp ){
			for (int ii=0;ii<nrows;ii++)
				for (int jj=0;jj<ncols;jj++)
					p[ii*ncols+jj]=U[ii*ncols+jj];			
				
			for (int s1=-1;s1<=1;s1+=2){
				for (int s2=-1;s2<=1;s2+=2){
					for (int i=(s1<0?nrows-1:0);(s1<0?i>=0:i<nrows);i+=s1){
						for (int j=(s2<0?ncols-1:0);(s2<0?j>=0:j<ncols);j+=s2){
							OneNode[0]=i;OneNode[1]=j;
							//if (notin(OneNode,src ,NumofSrc)){ //not a fixed node
							if (!(i==src[srcid*DIM+0] && j==src[srcid*DIM+1]) && (i+1)*dz<=Z[(ir*ncols+j)*3+0]){
								tmpU=eikslv2d(U,vp,OneNode,dx,dz,nrows,ncols);
								U[i*ncols+j]=MIN(U[i*ncols+j],tmpU);
							}
						}
					}
				}
			}
		}
			
		//sprintf(filename,"Src%d_downT2d.r%d",srcid,ir);writemod(U,filename,nrows,ncols);
		for (int i=0;i<nrows*ncols;i++) TTd[ir*nrows*ncols+i]=U[i];
		
		if (count>=MAXSwp) {printf("|\tPE_%d: ....Warning! downgoing FSM equation NOT Converged after %d sweeps for reflector %d...\n",pe,MAXSwp,ir+1);	
			sprintf(filename,"vpNotConverge4debug.tmp_PE%d",pe);writemod(vp,filename,nrows,ncols);
		}
		//End of downgoing...reinitiating...
		//if (pe==0) printf("*** %f ***Restart FSM\n",(double)(clock()-tmp)/CLOCKS_PER_SEC);
		tmp=clock();
		
		//inside model		
		for (int jj=1;jj<ncols-1;jj++){
			k2=-1;
			for (int kk=0;kk<nrows-1;kk++){
				if ((kk+1)*dz+dz<=Z[(ir*ncols+jj)*3+0] && (kk*dz+dz<=Z[(ir*ncols+jj-1)*3+0] || kk*dz+dz<=Z[(ir*ncols+jj+1)*3+0]))
					if (kk>=k2)
						k2=kk;
			}
				
			if (Z[(ir*ncols+jj)*3+0]<=nrows*dz && k2>0){ //find k2 and reflector is within the modeling area					
				if (k2*dz+dz<=Z[(ir*ncols+jj-1)*3+0] && k2*dz+dz<=Z[(ir*ncols+jj+1)*3+0])					
					{qx=(U[k2*ncols+jj+1]-U[k2*ncols+jj-1])/2.0/dx;
					if (U[k2*ncols+jj+1]>100 || U[k2*ncols+jj-1]>100) 
					printf("T[%d+1,%d]=%f, T[%d-1,%d]=%f\n",jj,k2,U[k2*ncols+jj+1],jj,k2,U[k2*ncols+jj-1]);
					}
				else if (k2*dz+dz<=Z[(ir*ncols+jj-1)*3+0])
					{qx=(U[k2*ncols+jj]-U[k2*ncols+jj-1])/dx;
					if (U[k2*ncols+jj]>100 || U[k2*ncols+jj-1]>100) 
					printf("T[%d,%d]=%f, T[%d-1,%d]=%f\n",jj,k2,U[k2*ncols+jj],jj,k2,U[k2*ncols+jj-1]);
					}
				else 
					{qx=(U[k2*ncols+jj+1]-U[k2*ncols+jj])/dx;
						if (U[k2*ncols+jj+1]>100 || U[k2*ncols+jj]>100) 
						printf("T[%d+1,%d]=%f, T[%d,%d]=%f\n",jj,k2,U[k2*ncols+jj+1],jj,k2,U[k2*ncols+jj]);
					}
									
				qz=(U[(k2+1)*ncols+jj]-U[k2*ncols+jj])/dz;
				dnx=Z[(ir*ncols+jj)*3+2];dnz=Z[(ir*ncols+jj)*3+1];
				rz=qz-2*(qx*dnx+qz*dnz)*dnz;
				U[(k2+1)*ncols+jj]=U[k2*ncols+jj]+qz*(Z[(ir*ncols+jj)*3+0]-k2*dz-dz)+rz*((k2+1)*dz+dz-Z[(ir*ncols+jj)*3+0]);
				k2vec[jj]=k2;
				if (U[(k2+1)*ncols+jj]<0) {
					printf("inside model: T[%d,%d]=%f\n",k2+1,jj,U[(k2+1)*ncols+jj]);
					printf("qx,qz,dnx,dnz,rz,T[%d],rd\n%f %f %f %f %f %f %f\n",k2,qx,qz,dnx,dnz,rz,U[k2*ncols+jj],Z[(ir*ncols+jj)*3+0]);exit(1);				
				}
			}	
			else if (Z[(ir*ncols+jj)*3+0]>nrows*dz)
				k2vec[jj]=nrows-1;
			else {printf("PE_%d: no k2 above reflector inside model...abort...\n",pe);exit(1);}
			
		}
		//if (pe==0) printf("*** %f ***Finish inside model\n",(double)(clock()-tmp)/CLOCKS_PER_SEC);
		tmp=clock();
		//at the 2 edges
		int jj=0;
		k2=-1;
		for (int kk=0;kk<nrows-1;kk++){
			if ((kk+1)*dz+dz<=Z[(ir*ncols+jj)*3+0] && kk*dz+dz<=Z[(ir*ncols+jj+1)*3+0])
				if (kk>=k2)
					k2=kk;
		}
			
		if (Z[(ir*ncols+jj)*3+0]<=nrows*dz && k2>0){ //find k2												
			qx=(U[k2*ncols+jj+1]-U[k2*ncols+jj])/dx;				
			qz=(U[(k2+1)*ncols+jj]-U[k2*ncols+jj])/dz;
			dnx=Z[(ir*ncols+jj)*3+2];dnz=Z[(ir*ncols+jj)*3+1];
			rz=qz-2*(qx*dnx+qz*dnz)*dnz;
			U[(k2+1)*ncols+jj]=U[k2*ncols+jj]+qz*(Z[(ir*ncols+jj)*3+0]-k2*dz-dz)+rz*((k2+1)*dz+dz-Z[(ir*ncols+jj)*3+0]);
			k2vec[jj]=k2;
			if (U[(k2+1)*ncols+jj]<0) {printf("at the left edge: T[%d,%d]=%f\n",k2+1,jj,U[(k2+1)*ncols+jj]);exit(1);}
		}
		else if (Z[(ir*ncols+jj)*3+0]>nrows*dz)
				k2vec[jj]=nrows-1;
		else { printf("PE_%d: no k2 above reflector at the edge ii=0...abort...\n",pe);exit(1);}				
		
		jj=ncols-1;		
		k2=-1;
		for (int kk=0;kk<nrows-1;kk++){
			if ((kk+1)*dz+dz<=Z[(ir*ncols+jj)*3+0] && kk*dz+dz<=Z[(ir*ncols+jj-1)*3+0])
				if (kk>=k2)
					k2=kk;
		}				
			
		if (Z[(ir*ncols+jj)*3+0]<=nrows*dz && k2>0){ //find k2								
			qx=(U[k2*ncols+jj]-U[k2*ncols+jj-1])/dx;				
			qz=(U[(k2+1)*ncols+jj]-U[k2*ncols+jj])/dz;
			dnx=Z[(ir*ncols+jj)*3+2];dnz=Z[(ir*ncols+jj)*3+1];
			rz=qz-2*(qx*dnx+qz*dnz)*dnz;
			U[(k2+1)*ncols+jj]=U[k2*ncols+jj]+qz*(Z[(ir*ncols+jj)*3+0]-k2*dz-dz)+rz*((k2+1)*dz+dz-Z[(ir*ncols+jj)*3+0]);
			k2vec[jj]=k2;
			if (U[(k2+1)*ncols+jj]<0) {printf("at the right edge: T[%d,%d]=%f\n",k2+1,jj,U[(k2+1)*ncols+jj]);exit(1);}
			//nodes except k2+1 reassiged to INF
			//for (int kk=0;kk<nz;kk++){
				//if (kk!=k2+1) U[kk*nrows*ncols+ii*ncols+jj]=INF;					}					
		}
		else if (Z[(ir*ncols+jj)*3+0]>nrows*dz)
				k2vec[jj]=nrows-1;
		else { printf("PE_%d: no k2 above reflector at the edge ii=nrows-1...abort...\n",pe);exit(1);}				
				
		//if (pe==0) {for (int jj=0;jj<ncols;jj++) printf("%f\n",(k2vec[jj]+1)*dz+dz-Z[(ir*ncols+jj)*3+0]);}
		//int tmp;cin>>tmp;
		/////////////////////////////
			
		tmp=clock();
		////////////////////
		for (int ii=0;ii<nrows;ii++){
			for (int jj=0;jj<ncols;jj++){
				if (ii<=k2vec[jj]) 
					U[ii*ncols+jj]=INF;
				p[ii*ncols+jj]=0.0;			
			}
		}
	
		//////////////////////
		//Start Sweepting
		count=0;
		while (!isconv(p,U,nrows,ncols) && count++<MAXSwp ){
			for (int ii=0;ii<nrows*ncols;ii++) p[ii]=U[ii];			
				
			for (int s1=-1;s1<=1;s1+=2){
				for (int s2=-1;s2<=1;s2+=2){
					for (int i=(s1<0?nrows-1:0);(s1<0?i>=0:i<nrows);i+=s1){
						for (int j=(s2<0?ncols-1:0);(s2<0?j>=0:j<ncols);j+=s2){
							OneNode[0]=i;OneNode[1]=j;
							//if (notin(OneNode,src ,NumofSrc)){ //not a fixed node
							if (i<=k2vec[j]) 
								{tmpU=eikslv2d(U,vp,OneNode,dx,dz,nrows,ncols);
								U[i*ncols+j]=MIN(U[i*ncols+j],tmpU);}							
							
						}
					}
				}
			}
		}
		
		//sprintf(filename,"Src%d_upT2d.r%d",srcid,ir);writemod(U,filename,nrows,ncols);
		for (int i=0;i<nrows*ncols;i++) TTu[ir*nrows*ncols+i]=U[i];
		
		if (count>=MAXSwp) printf("|\tPE_%d: ....Warning! upgoing FSM Equation NOT Converged after %d sweeps for reflector %d...\n",pe,MAXSwp,ir+1);	
		
		/*for (int i=0;i<srcid;i++) chr+=srnum[i]; 
		for (int idx=0;idx<srnum[srcid];idx++) {			
			int i=rec[(idx+chr)*DIM+0],j=rec[(idx+chr)*DIM+1];
			Tr[(ch0+idx)*(NumofR+1)+ir+1]=U[i*ncols+j];
		}*/
	}
	
	//if (pe==0) printf("*** %f ***Finish Reflection...Start Transmission...\n",(double)(clock()-tmp)/CLOCKS_PER_SEC);
	//////////////////////
	//Start Sweepting for transmission first arrivals
	for (int ii=0;ii<nrows;ii++){
		for (int jj=0;jj<ncols;jj++){
			if (ii==src[srcid*DIM+0] && jj==src[srcid*DIM+1]) 
				U[ii*ncols+jj]=0.0;
			else U[ii*ncols+jj]=INF;
			p[ii*ncols+jj]=0.0;			
		}
	}
	
	//////////////////////
	//Start Sweepting
	count=0;
	while (!isconv(p,U,nrows,ncols) && count++<MAXSwp ){
		for (int ii=0;ii<nrows;ii++)
			for (int jj=0;jj<ncols;jj++)
				p[ii*ncols+jj]=U[ii*ncols+jj];			
			
		for (int s1=-1;s1<=1;s1+=2){
			for (int s2=-1;s2<=1;s2+=2){
				for (int i=(s1<0?nrows-1:0);(s1<0?i>=0:i<nrows);i+=s1){
					for (int j=(s2<0?ncols-1:0);(s2<0?j>=0:j<ncols);j+=s2){
						OneNode[0]=i;OneNode[1]=j;
						//if (notin(OneNode,src ,NumofSrc)){ //not a fixed node
						if (!(i==src[srcid*DIM+0] && j==src[srcid*DIM+1])){
							tmpU=eikslv2d(U,vp,OneNode,dx,dz,nrows,ncols);
							U[i*ncols+j]=MIN(U[i*ncols+j],tmpU);
						}
					}
				}
			}
		}
	}
	/*for (int idx=0;idx<srnum[srcid];idx++) {
		int i=rec[(idx+chr)*DIM+0],j=rec[(idx+chr)*DIM+1];
		Tr[(ch0+idx)*(NumofR+1)+0]=U[i*ncols+j];
		//printf("output: ch0=%d,idx=%d,NumofR=%d, Tr=%f\n",ch0,idx,NumofR,Tr[(ch0+idx)*(NumofR+1)+0]);
	}*/
	//for (int idx=0;idx<srnum[srcid];idx++) {
		//printf("%f %f\n",Tr[(ch0+idx)*(NumofR+1)+0],Tr[(ch0+idx)*(NumofR+1)+1]);
	//}
	delete [] p;p=NULL;
	//if (count<MAXSwp)
		//printf("|\tPE_%d: Eikonal Equation Solver Converged after %d sweeps, %8.5f seconds\n\n",pe, count,(double)(clock()-tic)/CLOCKS_PER_SEC);
	//else printf("|\tPE_%d: ....Warning! Eikonal Equation Solver  NOT Converged after %d sweeps, %8.5f seconds...\n",pe,MAXSwp,(double)(clock()-tic)/CLOCKS_PER_SEC);
	if (count>=MAXSwp) printf("|\tPE_%d: ....Warning! Eikonal Equation Solver NOT Converged after %d sweeps, %8.5f seconds...\n",pe,MAXSwp,(double)(clock()-tic)/CLOCKS_PER_SEC);	
	
}

void usage()
{
	printf(" -----------------------------------------------------------------------------------------------------------------------\n");
	printf(" Wrong usage. See correct usage below:\n");
	printf(" mpirun -np number_of_processors RT_PFAST2D input_file_name [Maximum Memory Number=5] [Tomography Scheme=3, Joint] [>output file]\"\n");
	printf("      Inverse Modeling usage 1: \"mpirun -np 10 ../bin/RT_PFAST2D ./FSM2d01.inp\"\n");
	printf("      Inverse Modeling usage 2: \"mpirun -np 10 ../bin/RT_PFAST2D ./FSM2d01.inp 5 3 > FSM2D01.out\"\n");
	printf("      Forward Modeling usage 3: \"mpirun -np 10 ../bin/RT_PFAST2D_FW ./FSM2d01_FW.inp\"\n");
	printf("      Forward Modeling usage 4: \"mpirun -np 10 ../bin/RT_PFAST2D_FW_NoS ./FSM2d01_FW.inp > FSM2D01.out\"\n");
	printf(" -----------------------------------------------------------------------------------------------------------------------\n");

}

void info() 
{
	printf(" *******************************\n");
	printf(" * This is Travel Time Tomography program: Version 1.0            \n");
	printf(" * Parallel 2-D Traveltime Inversion based on Fast Sweeping Method and Adjoint-state Technique      \n");
	printf(" *                                                           \n");
	printf(" * written by  J.W. Huang                         \n");
	printf(" * See COPYING file for copying and redistribution conditions.\n");
	printf(" *******************************\n");
	printf(" \n");
}
#endif
