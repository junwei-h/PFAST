#ifndef _ADJOINT2D_METHODS
#define _ADJOINT2D_METHODS

#include "RT_pfsm2d.h"
#include "fftn.h"
#include <cmath>

#ifndef POS
	#define POS(x) (((x)+fabs(x))/2.0)
#endif

#ifndef NEG
	#define NEG(x) (((x)-fabs(x))/2.0) 
#endif

#ifndef MAXV
	#define MAXV 4500 
#endif

#ifndef MINV
	#define MINV 1600  
#endif

/*--------FFT or IFFT Switch----------*/
#define FORWARD_SCALE	0.0
#define INVERSE_SCALE	-1.0

using namespace std;

 void swap(double &a, double &b)
	{double dum=a; a=b; b=dum;}
	
int findin(int *node,int *CurrentNodes,int Num_CurrentNodes){
//to find the index of one node is in CurrentNodes
	//clock_t tic=clock();
	for (int i=0;i<Num_CurrentNodes;i++)
		if (CurrentNodes[i*DIM+0]==node[0] && CurrentNodes[i*DIM+1]==node[1] && CurrentNodes[i*DIM+2]==node[2])
			{//printf("findin function takes %8.5f seconds\n", (double)(clock()-tic)/CLOCKS_PER_SEC);
			return i;}
}

double adjoint_slv2d(double *U,double *T,int *OneNode,double dx,double dz, int nrows, int ncols){
	double a1,a2,b1,b2,La1,La2,Lb1,Lb2;
	
	if (T[OneNode[0]*ncols+OneNode[1]+1]>=INF || T[OneNode[0]*ncols+OneNode[1]]>=INF || T[OneNode[0]*ncols+OneNode[1]-1]>=INF \
	|| T[(OneNode[0]+1)*ncols+OneNode[1]]>=INF || T[(OneNode[0]-1)*ncols+OneNode[1]]>=INF)
		return 0.0;
	
	a1=POS(-(T[OneNode[0]*ncols+OneNode[1]+1]-T[OneNode[0]*ncols+OneNode[1]])/dx);
	a2=NEG(-(T[OneNode[0]*ncols+OneNode[1]]-T[OneNode[0]*ncols+OneNode[1]-1])/dx);
	b1=POS(-(T[(OneNode[0]+1)*ncols+OneNode[1]]-T[OneNode[0]*ncols+OneNode[1]])/dz);
	b2=NEG(-(T[OneNode[0]*ncols+OneNode[1]]-T[(OneNode[0]-1)*ncols+OneNode[1]])/dz);
	
	La1=POS(-(T[OneNode[0]*ncols+OneNode[1]]-T[OneNode[0]*ncols+OneNode[1]-1])/dx)*U[OneNode[0]*ncols+OneNode[1]-1];
	La2=NEG(-(T[OneNode[0]*ncols+OneNode[1]+1]-T[OneNode[0]*ncols+OneNode[1]])/dx)*U[OneNode[0]*ncols+OneNode[1]+1];
	Lb1=POS(-(T[OneNode[0]*ncols+OneNode[1]]-T[(OneNode[0]-1)*ncols+OneNode[1]])/dz)*U[(OneNode[0]-1)*ncols+OneNode[1]];
	Lb2=NEG(-(T[(OneNode[0]+1)*ncols+OneNode[1]]-T[OneNode[0]*ncols+OneNode[1]])/dz)*U[(OneNode[0]+1)*ncols+OneNode[1]];
	
	if (fabs((a1-a2)/dx+(b1-b2)/dz)>1e-9)
		return ((La1-La2)/dx+(Lb1-Lb2)/dz)/((a1-a2)/dx+(b1-b2)/dz);
	else{
		//printf("(%d,%d):a1=%f,a2=%f,b1=%f,b2=%f, (a1-a2+b1-b2)=%f\n",OneNode[0],OneNode[1],a1,a2,b1,b2,a1-a2+b1-b2);
		return 0.0;
	}
}

double Tadjoint_fsm2d(double *U,double *T,double *T_obs,double WR,int *rec, double *nxyz, int *src, int *srnum,int srcid,int MAXSwp,double dx,double dz,int nrows,int ncols,int pe){
	//printf("|\tPE_%d Called Adjoint Equation Solver for source #%d (%d,%d,%d)...\n",pe,srcid+1,src[srcid*DIM+0]+1,src[srcid*DIM+1]+1,src[srcid*DIM+2]+1);
	int count=0,OneNode[DIM],ch0=0,i,j;
	double nnx,nny,Tx,Tz,EE=0.0,Cr=1-WR;
	double *p=new double[nrows*ncols]; 
	int *rectag=new (nothrow)int[nrows*ncols]; //tag 1 if it is a receiver location, and 0 otherwise
	
	if (!p) {printf("PE_%d: Memory shortage in adjoint_fsm3d\n",pe);exit (1);}
	clock_t toc=clock();//start time recording
	////////////////////////////
	//Calculate Normal Vector at the surface defined by receiver locations
	//surfnorm(rec,NumofRec,nxyz,simds);
	
	///////////////////////////
	//Initialize lambda (U) boundary condition
	//clock_t timeini=clock();
	//printf("PE_%d: Starts Boundary Initialization...\n",pe);
	for (int k=0;k<nrows*ncols;k++){	U[k]=0.0;p[k]=0.0;	}
	for (int i=0;i<srcid;i++) ch0+=srnum[i];
	for (int k=0;k<nrows*ncols;k++) rectag[k]=0;
			
	for (int idx=0;idx<srnum[srcid];idx++) {
		if (T_obs[idx+ch0]<100) { //ONLY when transmission time picked
			i=rec[(idx+ch0)*DIM+0];j=rec[(idx+ch0)*DIM+1];
			rectag[i*ncols+j]=1;	
			EE+=Cr*(T_obs[idx+ch0]-T[i*ncols+j])*(T_obs[idx+ch0]-T[i*ncols+j])/2.0;
			if ( ! (rec[(idx+ch0)*DIM+0]==src[srcid*DIM+0] && rec[(idx+ch0)*DIM+1]==src[srcid*DIM+1]))// not at a source location
			{
				nny=nxyz[(idx+ch0)*DIM+0];nnx=nxyz[(idx+ch0)*DIM+1];
				//One sided difference for travel time at the boundary
					
				if (i+1>=nrows)
					Tz=(T[i*ncols+j]-T[(i-1)*ncols+j])/dz;
				else
					Tz=(T[(i+1)*ncols+j]-T[i*ncols+j])/dz;
						
				if (j+1>=ncols)
					Tx=(T[i*ncols+j]-T[i*ncols+j-1])/dx;
				else
					Tx=(T[i*ncols+j+1]-T[i*ncols+j])/dx;
							
				if (fabs(T_obs[idx+ch0]-T[i*ncols+j])>1e-9 && fabs(nnx*Tx+nny*Tz)>1e-6)
					U[i*ncols+j]=Cr*(T_obs[idx+ch0]-T[i*ncols+j])/(nnx*Tx+nny*Tz);
				else
					U[i*ncols+j]=0.0;
				//if (pe==0) printf("%d rectag(%d %d %d)=%d\n",idx+1,i+1,j+1,k+1,rectag[i*ncols+j]);
				//if (pe==0) printf("%d %d %d %d %15.9f\t%15.9f\t%15.9f\t%15.9f\t%15.9f\t%15.9f\t%15.9f\t%15.9f\t%15.9f\n",idx+1,i+1,j+1,k+1,U[i*ncols+j],T_obs[idx+ch0],T[i*ncols+j],nnx,nny,nnz,Tx,Tz,Tz);
			}
		}
	}
	//printf("PE_%d: %8.5f seconds for Boundary Initialization\n",pe,(double)(clock()-timeini)/CLOCKS_PER_SEC);
	/////////////////Done Initialization////////////////////
	 
	///////////////////////
	//Solve adjoint state	
	//while ( !ConvFlag && count++<10){
				
	while (count++<MAXSwp && !isconv(p,U,nrows,ncols)  ){
		//clock_t tic=clock();
		for (int kk=0;kk<nrows*ncols;kk++) p[kk]=U[kk];
		//printf("here\n");		
			for (int s1=-1;s1<=1;s1+=2){
				for (int s2=-1;s2<=1;s2+=2){
					//printf("PE_%d: Sweep (%d,%d,%d)...",pe,s3,s1,s2);clock_t tmp=clock();
					for ( i=(s1<0?nrows-2:1);(s1<0?i>=1:i<nrows-1);i+=s1){
						for ( j=(s2<0?ncols-2:1);(s2<0?j>=1:j<ncols-1);j+=s2){
							OneNode[0]=i;OneNode[1]=j;
							//it=search (vector_rec.begin(), vector_rec.end(), OneNode, OneNode+DIM);
							//if (notin(OneNode,rec,NumofRec)){ //not a fixed node 	
							if (fabs((double)rectag[i*ncols+j])<1e-9) {
								//printf("[s1 s2]=[%d %d]\t[i j]=[%d %d]\t[p q]=[%f %f]\n",s1,s2,i,j,p,q);
								//U[i*ncols+j]=q;//MAX(p,q);
								U[i*ncols+j]=adjoint_slv2d(U,T,OneNode,dx,dz,nrows,ncols);
							}
						}
					}
				
			//printf("takes %8.5f seconds\n", (double)(clock()-tmp)/CLOCKS_PER_SEC);
				}
			}
		//printf("Sweep #%d...takes time %8.5f seconds\n",count,(double)(clock()-tic)/CLOCKS_PER_SEC);	
	}
	delete [] p;p=NULL;
	delete [] rectag;rectag=NULL;
	//if (count<MAXSwp)
	//	printf("|\tPE_%d: Adjoint Equation Converged after %d sweeps, %8.5f seconds\n\n",pe,count,(double)(clock()-toc)/CLOCKS_PER_SEC);
	//else printf("|\tPE_%d: ....Warning! Adjoint Equation NOT Converged after %d sweeps, %8.5f seconds...\n",pe,MAXSwp,(double)(clock()-toc)/CLOCKS_PER_SEC);
	if (count>=MAXSwp) printf("|\tPE_%d: ....Warning!  Adjoint Equation of Transmission NOT Converged after %d sweeps, %8.5f seconds...\n",pe,MAXSwp,(double)(clock()-toc)/CLOCKS_PER_SEC);
	
	return EE; //time residual
}

double Radjoint_fsm2d(double *UR,double *Td, double *Tu,double *TR_obs,double WR,int *rec, double *nxyz, int *src, int *srnum, int srcid,int NumofR,int MAXSwp,double dx,double dz,int nrows,int ncols,int pe){
	//printf("|\tPE_%d Called Adjoint Equation Solver for source #%d (%d,%d,%d)...\n",pe,srcid+1,src[srcid*DIM+0]+1,src[srcid*DIM+1]+1,src[srcid*DIM+2]+1);
	int count=0,OneNode[DIM],ch0=0,i,j;
	double nnx,nny,Tx,Tz,EE=0.0,Cr=WR;
	double *p=new double[nrows*ncols],*U=new double[nrows*ncols]; 
	int *rectag=new (nothrow)int[nrows*ncols]; //tag 1 if it is a receiver location, and 0 otherwise
	//int *srctag=new (nothrow) int [nrows*ncols]; //tag 1 if it is a reflection point, and 0 otherwise.
	
	if (!p) {printf("PE_%d: Memory shortage in adjoint_fsm3d\n",pe);exit (1);}
	clock_t toc=clock();//start time recording
	
	///////////////////////////
	//Initialize lambda (U) boundary condition
	//clock_t timeini=clock();
	//printf("PE_%d: Starts Boundary Initialization...\n",pe);
	for (int k=0;k<nrows*ncols;k++){rectag[k]=0;UR[k]=0;}
	for (int i=0;i<srcid;i++) ch0+=srnum[i];
	
	for (int ir=0;ir<NumofR;ir++) {
		for (int k=0;k<nrows*ncols;k++) {U[k]=0.0;p[k]=0.0;}
		///////////////////
		//adjoint state of upgoing wave Tu
		for (int idx=0;idx<srnum[srcid];idx++) {
			if (TR_obs[(idx+ch0)*NumofR+ir]!=0) { //ONLY when reflection time picked 		
				i=rec[(idx+ch0)*DIM+0];j=rec[(idx+ch0)*DIM+1];rectag[i*ncols+j]=1;
			
				EE+=Cr*(TR_obs[(idx+ch0)*NumofR+ir]-Tu[ir*nrows*ncols+i*ncols+j])*(TR_obs[(idx+ch0)*NumofR+ir]-Tu[ir*nrows*ncols+i*ncols+j])/2.0;
				//if (k2vec[j]+1!=i)// receiver not overlap with a reflection point
				//{
					nny=nxyz[(idx+ch0)*DIM+0];nnx=nxyz[(idx+ch0)*DIM+1];
					//One sided difference for travel time at the boundary
						
					if (i+1>=nrows)
						Tz=(Tu[ir*nrows*ncols+i*ncols+j]-Tu[ir*nrows*ncols+(i-1)*ncols+j])/dz;
					else
						Tz=(Tu[ir*nrows*ncols+(i+1)*ncols+j]-Tu[ir*nrows*ncols+i*ncols+j])/dz;
							
					if (j+1>=ncols)
						Tx=(Tu[ir*nrows*ncols+i*ncols+j]-Tu[ir*nrows*ncols+i*ncols+j-1])/dx;
					else
						Tx=(Tu[ir*nrows*ncols+i*ncols+j+1]-Tu[ir*nrows*ncols+i*ncols+j])/dx;
								
					if (fabs(TR_obs[(idx+ch0)*NumofR+ir]-Tu[ir*nrows*ncols+i*ncols+j])>1e-9 && fabs(nnx*Tx+nny*Tz)>1e-6)
						U[i*ncols+j]=Cr*(TR_obs[(idx+ch0)*NumofR+ir]-Tu[ir*nrows*ncols+i*ncols+j])/(nnx*Tx+nny*Tz);
					else
						U[i*ncols+j]=0.0;
					//if (pe==0) printf("%d rectag(%d %d %d)=%d\n",idx+1,i+1,j+1,k+1,rectag[i*ncols+j]);
					//if (pe==0 && TR_obs[(idx+ch0)*NumofR+ir]!=-Tu[ir*nrows*ncols+i*ncols+j]) printf("%d %d %d %d %15.9f\t%15.9f\t%15.9f\t%15.9f\t%15.9f\t%15.9f\t%15.9f\n",idx+1,i+1,j+1,ir+1,U[i*ncols+j],TR_obs[(idx+ch0)*NumofR+ir],Tu[ir*nrows*ncols+i*ncols+j],nnx,nny,Tx,Tz);
				//}
			}
		}
		//printf("PE_%d: %8.5f seconds for Boundary Initialization\n",pe,(double)(clock()-timeini)/CLOCKS_PER_SEC);
		/////////////////Done Initialization////////////////////
		
		///////////////////////
		//Solve adjoint state	
		while (count++<MAXSwp && !isconv(p,U,nrows,ncols)  ){
			//clock_t tic=clock();
			for (int kk=0;kk<nrows*ncols;kk++) p[kk]=U[kk];
			//printf("here\n");		
				for (int s1=-1;s1<=1;s1+=2){
					for (int s2=-1;s2<=1;s2+=2){
						//printf("PE_%d: Sweep (%d,%d,%d)...",pe,s3,s1,s2);clock_t tmp=clock();
						for ( i=(s1<0?nrows-2:1);(s1<0?i>=1:i<nrows-1);i+=s1){
							for ( j=(s2<0?ncols-2:1);(s2<0?j>=1:j<ncols-1);j+=s2){
								OneNode[0]=i;OneNode[1]=j;
								//it=search (vector_rec.begin(), vector_rec.end(), OneNode, OneNode+DIM);
								//if (notin(OneNode,rec,NumofRec)){ //not a fixed node 	
								if (fabs((double)rectag[i*ncols+j])<1e-9) {
									//printf("[s1 s2]=[%d %d]\t[i j]=[%d %d]\t[p q]=[%f %f]\n",s1,s2,i,j,p,q);
									//U[i*ncols+j]=q;//MAX(p,q);
									U[i*ncols+j]=adjoint_slv2d(U,Tu+ir*nrows*ncols,OneNode,dx,dz,nrows,ncols);
								}
							}
						}
					
				//printf("takes %8.5f seconds\n", (double)(clock()-tmp)/CLOCKS_PER_SEC);
					}
				}
			//printf("Sweep #%d...takes time %8.5f seconds\n",count,(double)(clock()-tic)/CLOCKS_PER_SEC);	
		}
		if (count>=MAXSwp) printf("|\tPE_%d: ....Warning! Adjoint Equation of Reflected Wave NOT Converged after %d sweeps...\n",pe,MAXSwp);
		
		//for (int i=0;i<nrows*ncols;i++) UR[i]=UR[i]+U[i];	
		////////////////////////
		//adjoint state of downgoing Td
		//Solve adjoint state	
		/*for (int k=0;k<nrows*ncols;k++) {p[k]=0.0;}
		while (count++<MAXSwp && !isconv(p,U,nrows,ncols)  ){
			//clock_t tic=clock();
			for (int kk=0;kk<nrows*ncols;kk++) p[kk]=U[kk];
			//printf("here\n");		
				for (int s1=-1;s1<=1;s1+=2){
					for (int s2=-1;s2<=1;s2+=2){
						//printf("PE_%d: Sweep (%d,%d,%d)...",pe,s3,s1,s2);clock_t tmp=clock();
						for ( i=(s1<0?nrows-2:1);(s1<0?i>=1:i<nrows-1);i+=s1){
							for ( j=(s2<0?ncols-2:1);(s2<0?j>=1:j<ncols-1);j+=s2){
								OneNode[0]=i;OneNode[1]=j;
								//it=search (vector_rec.begin(), vector_rec.end(), OneNode, OneNode+DIM);
								//if (notin(OneNode,rec,NumofRec)){ //not a fixed node 	
								if (fabs((double)rectag[i*ncols+j])<1e-9) {
									//printf("[s1 s2]=[%d %d]\t[i j]=[%d %d]\t[p q]=[%f %f]\n",s1,s2,i,j,p,q);
									//U[i*ncols+j]=q;//MAX(p,q);
									U[i*ncols+j]=adjoint_slv2d(U,Td+ir*nrows*ncols,OneNode,dx,dz,nrows,ncols);
								}
							}
						}
					
				//printf("takes %8.5f seconds\n", (double)(clock()-tmp)/CLOCKS_PER_SEC);
					}
				}
			//printf("Sweep #%d...takes time %8.5f seconds\n",count,(double)(clock()-tic)/CLOCKS_PER_SEC);	
		}
		if (count>=MAXSwp) printf("|\tPE_%d: ....Warning! Adjoint Equation of Incident Wave NOT Converged after %d sweeps...\n",pe,MAXSwp);
		*/
		for (int i=0;i<nrows*ncols;i++) UR[i]=UR[i]+U[i];
	}
	delete [] p;p=NULL;
	delete [] rectag;rectag=NULL;
	delete [] U;U=NULL;
	//if (count<MAXSwp)
	//	printf("|\tPE_%d: Adjoint Equation Converged after %d sweeps, %8.5f seconds\n\n",pe,count,(double)(clock()-toc)/CLOCKS_PER_SEC);
	//else printf("|\tPE_%d: ....Warning! Adjoint Equation NOT Converged after %d sweeps, %8.5f seconds...\n",pe,MAXSwp,(double)(clock()-toc)/CLOCKS_PER_SEC);
	//if (count>=MAXSwp) printf("|\tPE_%d: ....Warning! Reflection Adjoint Equation NOT Converged after %d sweeps, %8.5f seconds...\n",pe,MAXSwp,(double)(clock()-toc)/CLOCKS_PER_SEC);
	
	return EE; //time residual
}

//Overload two functions: Tadjoint_fsm2d and Radjoint_fsm2d
double Tadjoint_fsm2d(double *U,double *T,double *T_obs,int *rec, double *nxyz, int *src, int *srnum,int srcid,int MAXSwp,double dx,double dz,int nrows,int ncols,int pe){
	//printf("|\tPE_%d Called Adjoint Equation Solver for source #%d (%d,%d,%d)...\n",pe,srcid+1,src[srcid*DIM+0]+1,src[srcid*DIM+1]+1,src[srcid*DIM+2]+1);
	//This is a overloaded function to handle spatial varying weighting factor
	
	int count=0,OneNode[DIM],ch0=0,i,j;
	double nnx,nny,Tx,Tz,EE=0.0,Cr=1.0;
	double *p=new double[nrows*ncols]; 
	int *rectag=new (nothrow)int[nrows*ncols]; //tag 1 if it is a receiver location, and 0 otherwise
	
	if (!p) {printf("PE_%d: Memory shortage in adjoint_fsm3d\n",pe);exit (1);}
	clock_t toc=clock();//start time recording
	////////////////////////////
	//Calculate Normal Vector at the surface defined by receiver locations
	//surfnorm(rec,NumofRec,nxyz,simds);
	
	///////////////////////////
	//Initialize lambda (U) boundary condition
	//clock_t timeini=clock();
	//printf("PE_%d: Starts Boundary Initialization...\n",pe);
	for (int k=0;k<nrows*ncols;k++){	U[k]=0.0;p[k]=0.0;	}
	for (int i=0;i<srcid;i++) ch0+=srnum[i];
	for (int k=0;k<nrows*ncols;k++) rectag[k]=0;
			
	for (int idx=0;idx<srnum[srcid];idx++) {
		if (T_obs[idx+ch0]<100) { //ONLY when transmission time picked
			i=rec[(idx+ch0)*DIM+0];j=rec[(idx+ch0)*DIM+1];
			rectag[i*ncols+j]=1;	
			EE+=Cr*(T_obs[idx+ch0]-T[i*ncols+j])*(T_obs[idx+ch0]-T[i*ncols+j])/2.0;
			if ( ! (rec[(idx+ch0)*DIM+0]==src[srcid*DIM+0] && rec[(idx+ch0)*DIM+1]==src[srcid*DIM+1]))// not at a source location
			{
				nny=nxyz[(idx+ch0)*DIM+0];nnx=nxyz[(idx+ch0)*DIM+1];
				//One sided difference for travel time at the boundary
					
				if (i+1>=nrows)
					Tz=(T[i*ncols+j]-T[(i-1)*ncols+j])/dz;
				else
					Tz=(T[(i+1)*ncols+j]-T[i*ncols+j])/dz;
						
				if (j+1>=ncols)
					Tx=(T[i*ncols+j]-T[i*ncols+j-1])/dx;
				else
					Tx=(T[i*ncols+j+1]-T[i*ncols+j])/dx;
							
				if (fabs(T_obs[idx+ch0]-T[i*ncols+j])>1e-9 && fabs(nnx*Tx+nny*Tz)>1e-6)
					U[i*ncols+j]=Cr*(T_obs[idx+ch0]-T[i*ncols+j])/(nnx*Tx+nny*Tz);
				else
					U[i*ncols+j]=0.0;
				//if (pe==0) printf("%d rectag(%d %d %d)=%d\n",idx+1,i+1,j+1,k+1,rectag[i*ncols+j]);
				//if (pe==0) printf("%d %d %d %d %15.9f\t%15.9f\t%15.9f\t%15.9f\t%15.9f\t%15.9f\t%15.9f\t%15.9f\t%15.9f\n",idx+1,i+1,j+1,k+1,U[i*ncols+j],T_obs[idx+ch0],T[i*ncols+j],nnx,nny,nnz,Tx,Tz,Tz);
			}
		}
	}
	//printf("PE_%d: %8.5f seconds for Boundary Initialization\n",pe,(double)(clock()-timeini)/CLOCKS_PER_SEC);
	/////////////////Done Initialization////////////////////
	 
	///////////////////////
	//Solve adjoint state	
	//while ( !ConvFlag && count++<10){
				
	while (count++<MAXSwp && !isconv(p,U,nrows,ncols)  ){
		//clock_t tic=clock();
		for (int kk=0;kk<nrows*ncols;kk++) p[kk]=U[kk];
		//printf("here\n");		
			for (int s1=-1;s1<=1;s1+=2){
				for (int s2=-1;s2<=1;s2+=2){
					//printf("PE_%d: Sweep (%d,%d,%d)...",pe,s3,s1,s2);clock_t tmp=clock();
					for ( i=(s1<0?nrows-2:1);(s1<0?i>=1:i<nrows-1);i+=s1){
						for ( j=(s2<0?ncols-2:1);(s2<0?j>=1:j<ncols-1);j+=s2){
							OneNode[0]=i;OneNode[1]=j;
							//it=search (vector_rec.begin(), vector_rec.end(), OneNode, OneNode+DIM);
							//if (notin(OneNode,rec,NumofRec)){ //not a fixed node 	
							if (fabs((double)rectag[i*ncols+j])<1e-9) {
								//printf("[s1 s2]=[%d %d]\t[i j]=[%d %d]\t[p q]=[%f %f]\n",s1,s2,i,j,p,q);
								//U[i*ncols+j]=q;//MAX(p,q);
								U[i*ncols+j]=adjoint_slv2d(U,T,OneNode,dx,dz,nrows,ncols);
							}
						}
					}
				
			//printf("takes %8.5f seconds\n", (double)(clock()-tmp)/CLOCKS_PER_SEC);
				}
			}
		//printf("Sweep #%d...takes time %8.5f seconds\n",count,(double)(clock()-tic)/CLOCKS_PER_SEC);	
	}
	delete [] p;p=NULL;
	delete [] rectag;rectag=NULL;
	//if (count<MAXSwp)
	//	printf("|\tPE_%d: Adjoint Equation Converged after %d sweeps, %8.5f seconds\n\n",pe,count,(double)(clock()-toc)/CLOCKS_PER_SEC);
	//else printf("|\tPE_%d: ....Warning! Adjoint Equation NOT Converged after %d sweeps, %8.5f seconds...\n",pe,MAXSwp,(double)(clock()-toc)/CLOCKS_PER_SEC);
	if (count>=MAXSwp) printf("|\tPE_%d: ....Warning!  Adjoint Equation of Transmission NOT Converged after %d sweeps, %8.5f seconds...\n",pe,MAXSwp,(double)(clock()-toc)/CLOCKS_PER_SEC);
	
	return EE; //time residual
}

double Radjoint_fsm2d(double *UR,double *Td, double *Tu,double *TR_obs,int *rec, double *nxyz, int *src, int *srnum, int srcid,int NumofR,int MAXSwp,double dx,double dz,int nrows,int ncols,int pe){
	//printf("|\tPE_%d Called Adjoint Equation Solver for source #%d (%d,%d,%d)...\n",pe,srcid+1,src[srcid*DIM+0]+1,src[srcid*DIM+1]+1,src[srcid*DIM+2]+1);
	//This is a overloaded function to handle spatial varying weighting factor
	int count=0,OneNode[DIM],ch0=0,i,j;
	double nnx,nny,Tx,Tz,EE=0.0,Cr=1.0;
	double *p=new double[nrows*ncols],*U=new double[nrows*ncols]; 
	int *rectag=new (nothrow)int[nrows*ncols]; //tag 1 if it is a receiver location, and 0 otherwise
	//int *srctag=new (nothrow) int [nrows*ncols]; //tag 1 if it is a reflection point, and 0 otherwise.
	
	if (!p) {printf("PE_%d: Memory shortage in adjoint_fsm3d\n",pe);exit (1);}
	clock_t toc=clock();//start time recording
	
	///////////////////////////
	//Initialize lambda (U) boundary condition
	//clock_t timeini=clock();
	//printf("PE_%d: Starts Boundary Initialization...\n",pe);
	for (int k=0;k<nrows*ncols;k++){rectag[k]=0;UR[k]=0;}
	for (int i=0;i<srcid;i++) ch0+=srnum[i];
	
	for (int ir=0;ir<NumofR;ir++) {
		for (int k=0;k<nrows*ncols;k++) {U[k]=0.0;p[k]=0.0;}
		///////////////////
		//adjoint state of upgoing wave Tu
		for (int idx=0;idx<srnum[srcid];idx++) {
			if (TR_obs[(idx+ch0)*NumofR+ir]!=0) { //ONLY when reflection time picked 		
				i=rec[(idx+ch0)*DIM+0];j=rec[(idx+ch0)*DIM+1];rectag[i*ncols+j]=1;
			
				EE+=Cr*(TR_obs[(idx+ch0)*NumofR+ir]-Tu[ir*nrows*ncols+i*ncols+j])*(TR_obs[(idx+ch0)*NumofR+ir]-Tu[ir*nrows*ncols+i*ncols+j])/2.0;
				//if (k2vec[j]+1!=i)// receiver not overlap with a reflection point
				//{
					nny=nxyz[(idx+ch0)*DIM+0];nnx=nxyz[(idx+ch0)*DIM+1];
					//One sided difference for travel time at the boundary
						
					if (i+1>=nrows)
						Tz=(Tu[ir*nrows*ncols+i*ncols+j]-Tu[ir*nrows*ncols+(i-1)*ncols+j])/dz;
					else
						Tz=(Tu[ir*nrows*ncols+(i+1)*ncols+j]-Tu[ir*nrows*ncols+i*ncols+j])/dz;
							
					if (j+1>=ncols)
						Tx=(Tu[ir*nrows*ncols+i*ncols+j]-Tu[ir*nrows*ncols+i*ncols+j-1])/dx;
					else
						Tx=(Tu[ir*nrows*ncols+i*ncols+j+1]-Tu[ir*nrows*ncols+i*ncols+j])/dx;
								
					if (fabs(TR_obs[(idx+ch0)*NumofR+ir]-Tu[ir*nrows*ncols+i*ncols+j])>1e-9 && fabs(nnx*Tx+nny*Tz)>1e-6)
						U[i*ncols+j]=Cr*(TR_obs[(idx+ch0)*NumofR+ir]-Tu[ir*nrows*ncols+i*ncols+j])/(nnx*Tx+nny*Tz);
					else
						U[i*ncols+j]=0.0;
					//if (pe==0) printf("%d rectag(%d %d %d)=%d\n",idx+1,i+1,j+1,k+1,rectag[i*ncols+j]);
					//if (pe==0 && TR_obs[(idx+ch0)*NumofR+ir]!=-Tu[ir*nrows*ncols+i*ncols+j]) printf("%d %d %d %d %15.9f\t%15.9f\t%15.9f\t%15.9f\t%15.9f\t%15.9f\t%15.9f\n",idx+1,i+1,j+1,ir+1,U[i*ncols+j],TR_obs[(idx+ch0)*NumofR+ir],Tu[ir*nrows*ncols+i*ncols+j],nnx,nny,Tx,Tz);
				//}
			}
		}
		//printf("PE_%d: %8.5f seconds for Boundary Initialization\n",pe,(double)(clock()-timeini)/CLOCKS_PER_SEC);
		/////////////////Done Initialization////////////////////
		
		///////////////////////
		//Solve adjoint state	
		while (count++<MAXSwp && !isconv(p,U,nrows,ncols)  ){
			//clock_t tic=clock();
			for (int kk=0;kk<nrows*ncols;kk++) p[kk]=U[kk];
			//printf("here\n");		
				for (int s1=-1;s1<=1;s1+=2){
					for (int s2=-1;s2<=1;s2+=2){
						//printf("PE_%d: Sweep (%d,%d,%d)...",pe,s3,s1,s2);clock_t tmp=clock();
						for ( i=(s1<0?nrows-2:1);(s1<0?i>=1:i<nrows-1);i+=s1){
							for ( j=(s2<0?ncols-2:1);(s2<0?j>=1:j<ncols-1);j+=s2){
								OneNode[0]=i;OneNode[1]=j;
								//it=search (vector_rec.begin(), vector_rec.end(), OneNode, OneNode+DIM);
								//if (notin(OneNode,rec,NumofRec)){ //not a fixed node 	
								if (fabs((double)rectag[i*ncols+j])<1e-9) {
									//printf("[s1 s2]=[%d %d]\t[i j]=[%d %d]\t[p q]=[%f %f]\n",s1,s2,i,j,p,q);
									//U[i*ncols+j]=q;//MAX(p,q);
									U[i*ncols+j]=adjoint_slv2d(U,Tu+ir*nrows*ncols,OneNode,dx,dz,nrows,ncols);
								}
							}
						}
					
				//printf("takes %8.5f seconds\n", (double)(clock()-tmp)/CLOCKS_PER_SEC);
					}
				}
			//printf("Sweep #%d...takes time %8.5f seconds\n",count,(double)(clock()-tic)/CLOCKS_PER_SEC);	
		}
		if (count>=MAXSwp) printf("|\tPE_%d: ....Warning! Adjoint Equation of Reflected Wave NOT Converged after %d sweeps...\n",pe,MAXSwp);
		
		for (int i=0;i<nrows*ncols;i++) UR[i]=UR[i]+U[i];
	}
	delete [] p;p=NULL;
	delete [] rectag;rectag=NULL;
	delete [] U;U=NULL;
	
	return EE; //time residual
}

bool check_conv(double *dvp,double stepsize,double *vp,double *dEE, int Iterk, double dt, int nrows, int ncols){
	double mynorm=0,chgE;
	
	if (Iterk>=2) chgE=(dEE[Iterk-2]-dEE[Iterk-1])/dEE[Iterk-2];
	else chgE=1;
	
	for (int k=0;k<nrows*ncols;k++){
		mynorm=(stepsize*stepsize*dvp[k]*dvp[k]/vp[k]/vp[k]>=mynorm)?\
		(stepsize*stepsize*dvp[k]*dvp[k]/vp[k]/vp[k]):mynorm;
				//mynorm=mynorm+dm[k]*dm[k]/vp[k]*vp[k];
	}
	mynorm=sqrt(mynorm);
	//if (mynorm<=eps || sqrt(2*dEE[Iterk-1])*1000<=dt || chgE <=0.01) //Total update amount < eps or residual error less than dt ms
	if (chgE <=0.001)	
		return true;
	else if ((mynorm/nrows/ncols)>=1) // average update is as twice as the previous velocity is considered as over corrected
		{printf("Over corrected..abort...\n");exit(1);}
	else
		return false;
}

void updatevp(double *vp,double *dm, double stepsize, int nrows, int ncols){
	for (int k=0;k<nrows*ncols;k++){
		if (fabs(vp[k]-330)>1e-9) {
			vp[k]=vp[k]+stepsize*dm[k];
			if (vp[k]>MAXV) vp[k]=MAXV;
			if (vp[k]<MINV)	vp[k]=MINV;
		}
	}
}

void regu2d(double *Dj,double *vp,double nux, double nuz, double dx,double dz, int nrows, int ncols){
	int ret,dims2[]={ncols,nrows};
	double *rDj=new double [nrows*ncols],*iDj=new double [nrows*ncols];
	double *kx=new double [ncols],*kz=new double [nrows];
	double kxmax,kzmax,kxs,kzs;
	
	for (int i=0;i<nrows;i++) {
		for (int j=0;j<ncols;j++){
			rDj[i*ncols+j]=Dj[i*ncols+j];iDj[i*ncols+j]=0.0;
		}
	}
	ret = fftn(2, dims2, rDj, iDj,  -1, FORWARD_SCALE); //combination of (-1, 0.0)==fftn.m in matlab
	if (ret) exit(1);
	
	kxmax=PI/dx;kzmax=PI/dz;kxs=2*PI/ncols/dx;kzs=2*PI/nrows/dz;
	
	for(int j=0;j<ncols;j++)
		kx[j]=-kxmax+kxs*j;
	for(int i=0;i<nrows;i++)	
		kz[i]=-kzmax+kzs*i;
	
	for(int j=0;j<ncols/2;j++)
		swap(kx[j],kx[ncols/2+j]);
	for(int i=0;i<nrows/2;i++)	
		swap(kz[i],kz[nrows/2+i]);
		
	for (int i=0;i<nrows;i++){
		for (int j=0;j<ncols;j++){
			rDj[i*ncols+j]=rDj[i*ncols+j]/(1+nux*nux*kx[j]*kx[j]+nuz*nuz*kz[i]*kz[i]);
			iDj[i*ncols+j]=iDj[i*ncols+j]/(1+nux*nux*kx[j]*kx[j]+nuz*nuz*kz[i]*kz[i]);		
		}
	}
	ret = fftn(2, dims2, rDj, iDj,  1, INVERSE_SCALE); //combination of (1, -1.0)==ifftn.m in matlab
	if (ret) exit(1);
	
	for (int i=0;i<nrows;i++) {
		for (int j=0;j<ncols;j++) {
			//if (fabs(vp[i*ncols+j]-330)<1e-9 || i==0 || i==nrows-1 || j==0 || j==ncols-1) 
			if (fabs(vp[i*ncols+j]-330)<1e-9 )
			//if (fabs(rDj[i*ncols+j])<1e-10)
				Dj[i*ncols+j]=0.0;
			else
				Dj[i*ncols+j]=rDj[i*ncols+j];
		}
	}
	delete [] rDj;rDj=NULL;
	delete [] iDj;iDj=NULL;
	delete [] kx;kx=NULL;
	delete [] kz;kz=NULL;	
}

double beta_HS(double *gk1,double *gk0,double *dk, int nrows, int ncols){
	double beta=0.0,tmpa=0.0,tmpb=0.0;
	
	for (int i=0;i<nrows;i++){
		for (int j=0;j<ncols;j++){
			tmpa+=gk1[i*ncols+j]*(gk1[i*ncols+j]-gk0[i*ncols+j]);
			tmpb+=-dk[i*ncols+j]*(gk1[i*ncols+j]-gk0[i*ncols+j]);
		}
	}
	beta=tmpa/tmpb;
	
	return beta;
}

double beta_CGDESCENT(double *gk1,double *gk0,double *dk, int nrows, int ncols){
	double beta=0.0,dTy=0.0,ny=0.0;
	
	for (int i=0;i<nrows*ncols;i++){
		dTy+=dk[i]*(-gk1[i]+gk0[i]);
		ny+=(-gk1[i]+gk0[i])*(-gk1[i]+gk0[i]);
	}
	for (int i=0;i<nrows*ncols;i++){
		beta+=((-gk1[i]+gk0[i])-2*dk[i]*ny/dTy)*(-gk1[i])/dTy;
	}
	return beta;
}

double gd2sk(double *gk1,double *gk0,double *dk, double alpha1, double *rhok, double *sk, double *yk, int Iterk, int lm, int nrows, int ncols){
	double gamma=0.0,dTy=0.0,ny=0.0;
	
	for (int i=0;i<nrows*ncols;i++){
		dTy+=alpha1*dk[i]*(-gk1[i]+gk0[i]);
		ny+=(-gk1[i]+gk0[i])*(-gk1[i]+gk0[i]);
	}
	gamma=dTy/ny;
	rhok[(Iterk-2)%lm]=1.0/dTy;
	for (int i=0;i<nrows*ncols;i++){
		sk[((Iterk-2)%lm)*nrows*ncols+i]=alpha1*dk[i];
		yk[((Iterk-2)%lm)*nrows*ncols+i]=-gk1[i]+gk0[i];
	}
	return gamma;
}
/*
void RTweight(double *wT, double *lambdaRT,double *lambdaR,double *lambdaT,double *loc_E,double *loc_Er,double *loc_Et,int nrows,int ncols)
{
	for (int kk=0;kk<nrows*ncols;kk++) //lambdaRT[kk]=wT*lambdaR[kk]+lambdaT[kk];
		lambdaRT[kk]=wT[kk]*lambdaR[kk]+lambdaT[kk];
			//loc_E[0]=wT*loc_Er[0]+loc_Et[0];
			loc_E[0]=loc_Et[0]+loc_Er[0];
}

double zoom(double alpha0,double alpha1,double phi0,double phi1,double dphi0,double dphi1){
	double a,b,alphac;
	
	a=(2*(phi0-phi1)+(alpha1-alpha0)*(dphi0-dphi1))/((alpha1-alpha0)*(alpha1-alpha0)*(alpha1-alpha0));
	b=(3*(phi1-phi0)-(alpha1-alpha0)*(dphi1+2*dphi0))/((alpha1-alpha0)*(alpha1-alpha0));
	
	if ((b*b-3*a*dphi0)>=0) alphac=alpha0+(-b+sqrt(b*b-3*a*dphi0))/3.0/a;
	else alphac=alpha1-(alpha1-alpha0)*dphi1/(dphi1-dphi0);
	return alphac;
}
*/
double cubicInterp(double alpha0,double alpha1,double phi0,double phi1,double dphi0,double dphi1){
	double a,b,c,d;
	a=(-2.*(phi1-phi0)+(dphi0+dphi1)*(alpha1-alpha0))/((alpha1-alpha0)*(alpha1-alpha0)*(alpha1-alpha0));
	b=(3.*(phi1-phi0)-(2.*dphi0+dphi1)*(alpha1-alpha0))/((alpha1-alpha0)*(alpha1-alpha0));
	c=dphi0;d=phi0;
	
	if (fabs(a)>0 && (b*b-3*a*c)>0)
		return alpha0+(-b+sqrt(b*b-3*a*c))/(3.*a);
	if (a==0 && fabs(b)>0 && (b*b-3*a*c)<0)
		return -c/(2.*b);
}

double magnitudeOf(double *g,int nrows,int ncols) {
	double maxg=abs(g[0]);
	for (int kk=1;kk<nrows*ncols;kk++) 
		if (abs(g[kk])>maxg) maxg=abs(g[kk]);
	
	return 1./maxg;
}

double meanOf(double *g,int nrows,int ncols) {
	double meang=0;
	for (int kk=0;kk<nrows*ncols;kk++) 
		meang=meang+g[kk];
	meang=meang/(nrows*ncols);
	
	return 1./fabs(meang);
}

/*double cubicInterp(double alpha0,double alpha1,double phi0,double phi1,double dphi0,double dphi1){
	MatDoub AMat(4,4);
	VecDoub XVec(4),BVec(4);
	double a1,a2,atmp;
	AMat[0][0]=alpha0*alpha0*alpha0;AMat[0][1]=alpha0*alpha0;AMat[0][2]=alpha0;AMat[0][3]=1.;
	AMat[1][0]=alpha1*alpha1*alpha1;AMat[1][1]=alpha1*alpha1;AMat[1][2]=alpha1;AMat[1][3]=1.;
	AMat[2][0]=3.*alpha0*alpha0;AMat[2][1]=2.*alpha0;AMat[2][2]=1.;AMat[2][3]=0.;
	AMat[3][0]=3.*alpha1*alpha1;AMat[3][1]=2.*alpha1;AMat[3][2]=1.;AMat[3][3]=0.;
	
	BVec[0]=phi0;BVec[1]=phi1;BVec[2]=dphi0;BVec[3]=dphi1;
	
	LUdcmp alu(AMat);
	alu.solve(BVec,XVec);
	//printf("XVec=[%g %g %g %g]\n",XVec[0],XVec[1],XVec[2],XVec[3]);
	
	if ((XVec[1]*XVec[1]-4.*XVec[0]*XVec[2])>=0) {
		a1=(-XVec[1]+sqrt(XVec[1]*XVec[1]-4.*XVec[0]*XVec[2]))/(3.*XVec[0]);
		a2=(-XVec[1]-sqrt(XVec[1]*XVec[1]-4.*XVec[0]*XVec[2]))/(3.*XVec[0]);	
		if (a1<0 && a2<0) {
			printf("!!!Warning...Cubic Interpolation Found No Posotive Step, One Iteration of Secant Method is Used Instead!!!\n");
			return alpha1-(alpha1-alpha0)*dphi1/(dphi1-dphi0);
		}
		else if ((6.*XVec[0]*MAX(a1,a2)+2.*XVec[1])>0)
			return MAX(a1,a2);		
	}
	else {
		printf("!!!Warning...Cubic Interpolation Failed, One Iteration of Secant Method is Used Instead!!!\n");
		return alpha1-(alpha1-alpha0)*dphi1/(dphi1-dphi0);
	}
}*/

void RTinfo(int RTtag){
	switch(RTtag){
		case 3:
			printf("====================================================================\n");
			printf("||This is a 2D Joint Inversion of Reflection and Transmission Time!||\n");
			printf("====================================================================\n");
		break;
				
		case 1:
			printf("=================================================================\n");
			printf("||        This is a 2D Inversion of Transmission Time ONLY!      ||\n");
			printf("=================================================================\n");
		break;
				
		case 2:
			printf("=================================================================\n");
			printf("||		This is a 2D Inversion of Reflection Time ONLY!		   ||\n");
			printf("=================================================================\n");
		break;
	}
}
#endif
