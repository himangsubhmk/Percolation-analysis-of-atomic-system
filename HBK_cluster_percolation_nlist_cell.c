/**
   @Program : cluster.c
   @Note    : Find cluster of active atom for plastic event
   @Note    : Atom less than 2.06A (first minima of gr) => same cluster
 **/
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>

#define totframe 1000
#define del_rc 0.1




int main(int argc, char *argv[])
{     

    FILE *fp,*fp1;
    char FILE[512],FILE1[512];

    double L;
    int N;  
    long i,j,k,l;
    double t;

    double charge,x_cor,y_cor,z_cor,tild,act,max_act;
    double xr,yr,zr;
    int iframe,id,itype,nx,ny,nz;
    char char_dummy[100];
    double dummy_double;

    int index;  
    double Li,Lf,dr;

    double gmax;
    int sample, cycle;
    int timestep;

    double **xyz_0,**xyz_1;
    double *activity, *max_activity;
    char *accum_string, *Dcalc_string;
    double cutoff;

    double **xsmax;
    int smax,smaxid;
    
    //------------------------SM update----------------------------------
  
    gmax=0.069;
    sample=1;
    cycle=atoi(argv[1]);
    cutoff=atof(argv[2]);
    accum_string=argv[3]; // options: accumulated, not_accumulated
    Dcalc_string=argv[4]; // options: Dcyc, Dinit
  
    if (strcmp(Dcalc_string,"Dinit")==0){
        sprintf(FILE, "/Data4/swarnendu/FatigueFailure/analysis/damage/finiterate/N64K/EIS-7.00/D2min_from_initial_config/make_dump/gmax0.069/gmax_g0.069_s1_c%d.dat",cycle);
    }else if (strcmp(Dcalc_string,"Dcyc")==0){
      //sprintf(FILE, "/Data4/swarnendu/FatigueFailure/properties/D2min/finiterate/N64K/EIS-7.00/make_dump/gmax0.069/gmax_g0.069_s1_c%d.dat",cycle);
      sprintf(FILE,"problamatic_conf_EIS-6.96/gmax_g0.069_s1_c%d.dat",cycle);
    }else{
        printf("4th argument should be Dcyc or Dinit\n");
        return 0;
    }
    fp=fopen(FILE,"r");
    while(!feof(fp)){
        fscanf(fp, "%*[^\n]\n");
        fscanf(fp,"%d\n",&timestep);
        fscanf(fp, "%*[^\n]\n");
        fscanf(fp,"%d\n",&N);

        xyz_0=malloc(N*sizeof(double *));
	xsmax=malloc(N*sizeof(double *));
        xyz_1=malloc(N*sizeof(double *));
        activity=malloc(N*sizeof(double));
        max_activity=malloc(N*sizeof(double));
        for(j=0;j<N;j++){
            xyz_0[j]=malloc(3*sizeof(double));
            xyz_1[j]=malloc(3*sizeof(double));
	    xsmax[j]=malloc(3*sizeof(double));
        }
        for(j=0;j<N;j++){
            activity[j]=0.0;
            max_activity[j]=0.0;
            for(k=0;k<3;k++){
                xyz_0[j][k]=0.0;
                xyz_1[j][k]=0.0;
		xsmax[j][k]=0.0;
            }
        }
        fscanf(fp, "%*[^\n]\n");
        fscanf(fp,"%lf %lf \n",&Li,&Lf);
        L=Lf-Li;
        fscanf(fp, "%*[^\n]\n");
        fscanf(fp, "%*[^\n]\n");
        fscanf(fp, "%*[^\n]\n");
        fscanf(fp, "%*[^\n]\n");
        for (j=0;j<N;j++){
            fscanf(fp, "%lf %lf %lf %lf %lf %*f\n", &x_cor,&y_cor,&z_cor,&act,&max_act);
        	xyz_0[j][0]=x_cor;
	        xyz_0[j][1]=y_cor;
	        xyz_0[j][2]=z_cor;
            activity[j]=act;
            max_activity[j]=max_act;
        }
    }
    fclose(fp);
    tild=0;

   printf("File read completed\n");
   

    int iplist,*PLIST,*LIST,*ID,*ns,*tot_ns;
    PLIST=malloc(N*sizeof(int));
    LIST =malloc(N*sizeof(int));
    ID   =malloc(N*sizeof(int));
    ns   =malloc((N+1)*sizeof(int));
    tot_ns   =malloc((N+1)*sizeof(int));
    int iroot,inew,n,ncl,k1,k2;
    
    for(i=0;i<=N;i++)tot_ns[i]=0;

    iplist=0;
    if (strcmp(accum_string,"not_accumulated")==0){
        for (i=0;i<N;i++){
            if (activity[i]>cutoff) PLIST[iplist++]=i;
        }
    }else if (strcmp(accum_string,"accumulated")==0){
         for (i=0;i<N;i++){
            if (max_activity[i]>cutoff) PLIST[iplist++]=i;
        } 
    }else{
        printf("3rd argument should be accumulated or not_accumulated\n");
        return 0;
    }

    
      /* //==========Transform co-ordinate to wrapped one========== */
      for(i=0;i<N;i++){
       for(j=0;j<3;j++){
	while(xyz_0[i][j]<0)xyz_0[i][j]+=L;
	while(xyz_0[i][j]>L)xyz_0[i][j]-=L;
       }
      }
      /* //========================================================= */

    printf("total_active=%d\n",iplist);
    if(iplist==0){
        printf("No active particles\n");
        return 0;
    }

    int ntot=0;
    ncl=0;
    smax=0;


    for(i=0;i<N;i++)ID[i]=-1;
    for(i=0;i<=N;i++)ns[i]=0;


    for(i=0;i<iplist;i++){

        iroot=PLIST[i];

        if(ID[iroot]==-1){
          ncl++;
          ID[iroot]=ncl;
          n=1;
          LIST[n]=iroot;
          k1=n;
          k2=k1;
          while(k1<=k2){
          for(k=k1;k<=k2;k++){
            iroot=LIST[k];
            for(j=0;j<iplist;j++){
              inew=PLIST[j];
              
              if(ID[inew]==-1){

            xr=xyz_0[iroot][0]-xyz_0[inew][0];
            yr=xyz_0[iroot][1]-xyz_0[inew][1];
            zr=xyz_0[iroot][2]-xyz_0[inew][2];
            
            //xr=xr-tild*round(zr/L);
            xr=xr-L*round(xr/L);
            yr=yr-L*round(yr/L);
            zr=zr-L*round(zr/L);
            dr=sqrt(xr*xr+yr*yr+zr*zr);
            
            
            if(dr<1.4){
              ID[inew]=ncl;
              n++;
              LIST[n]=inew;
            }
              }
            }//loop j
            
          }//loop k

            
          k1=k2+1;
          k2=n;
          }//while

	  
	  if(n>smax){
	    smax=n;
	    smaxid=ncl;
	  }

	  
          ns[n]++;
          //printf("ncl=%d n=%d\n",ncl,n);
          ntot+=n;
          
        }//if
	
    }//loop i

    l=0;
    for(i=0;i<iplist;i++){
      j=PLIST[i];
      if(ID[j]==smaxid){
	for(k=0;k<3;k++)
	  xsmax[l][k]=xyz_0[j][k];
	l++;
      }
    }
    //printf("l=%d smax=%d\n",l,smax);


    int Pdir;
    
    int percolation=1;
    if(percolation){

    
    
     //for image cluster counting
    int ncopy,ncopy3;
    ncopy=3;
    ncopy3=ncopy*ncopy*ncopy;
    double **img_xyz_0;
    img_xyz_0=malloc(ncopy3*smax*sizeof(double *));
    for(j=0;j<ncopy3*smax;j++)img_xyz_0[j]=malloc(3*sizeof(double));

    int img_iplist,*img_LIST,*img_ID,*img_ns,*img_tot_ns;

    img_LIST =malloc(ncopy3*smax*sizeof(int));
    img_ID   =malloc(ncopy3*smax*sizeof(int));
    img_ns   =malloc((ncopy3*smax+1)*sizeof(int));
    img_tot_ns   =malloc((ncopy3*smax+1)*sizeof(int));
    int img_smax;
    int  Ncopy;
    double temp,Lcopy,iLcopy,dgamma,tildcopy;
    double P,Pinf,PP,Pbin;


    
    
      // Heading to cluster counting considering image
      
      dgamma=tild/L;
      Ncopy=ncopy3*smax;
      Lcopy=ncopy*L;
      iLcopy=1.0/Lcopy;
      tildcopy=ncopy*tild;
      img_smax=0.0;
      img_iplist=0;

      index=0;
      for(i=0;i<smax;i++){

      for(nx=0;nx<ncopy;nx++)
      	for(ny=0;ny<ncopy;ny++)
      	  for(nz=0;nz<ncopy;nz++){
      	    temp=xsmax[i][0]-dgamma*xsmax[i][2];
      	    img_xyz_0[index][0]=temp+nx*L;
      	    img_xyz_0[index][1]=xsmax[i][1]+ny*L;
      	    img_xyz_0[index][2]=xsmax[i][2]+nz*L;

      	    img_xyz_0[index][0]+=dgamma*img_xyz_0[index][2];
	    
	    index++;
      	  }
      }
      printf("img_iplist=%d\n",img_iplist);
      //if(img_iplist==0)continue;

      
      //========== Get the neighbour list ==========================/
      //18Aug24: Doing neighbourlist for large N has N^2 computation;
      //Here I try to reduce that by making cell/calling particles
      //from that cell and neighbouring cell only to do the neighbour list...
      //lets see how it works...
      
      //19Aug24: It works very well; If you use this for other code
      //be careful about cutoff 1.5 and 1.4, also number density

      double drcut=1.4; //distancecutoff for considering neighbour
      double rho=1.2; //number density
      double rho2=2.0*rho;
      int maxneighbour=rho2*4/3*M_PI*drcut*drcut*drcut;
      int **neighbourLIST;
      neighbourLIST=malloc(Ncopy*sizeof(int *));
      for(j=0;j<Ncopy;j++)neighbourLIST[j]=malloc(maxneighbour*sizeof(int));
      int *nlist;
      nlist=malloc(Ncopy*sizeof(int));
      for(i=0;i< Ncopy;i++)nlist[i]=0;
      
      double approxdecllwidth=drcut+0.1;//drcut is the magic distance,so...
      
      int maxcell=(int)(Lcopy/approxdecllwidth);
      printf("maxcell=%d\n",maxcell);
      int maxcell3=maxcell*maxcell*maxcell;
      double cellwidth=Lcopy/maxcell;
      printf("cellwidth=%e\n",cellwidth);
      int maxcellpart=(rho2*pow(cellwidth,3)); //2*rho for safety 
      printf("maxcellpart=%d\n",maxcellpart);
      int **CELLLIST,*clist;
      CELLLIST=malloc(maxcell3*sizeof(int *));
      for(j=0;j<maxcell3;j++)
	  CELLLIST[j]=malloc(maxcellpart*sizeof(int));
      clist=malloc(maxcell3*sizeof(int));
      for(i=0;i<maxcell3;i++)clist[i]=0;
      

      
      for(int ipart=0;ipart< Ncopy;ipart++){
	int ix=img_xyz_0[ipart][0]/cellwidth;
	int iy=img_xyz_0[ipart][1]/cellwidth;
	int iz=img_xyz_0[ipart][2]/cellwidth;

	//printf("ix %d iy=%d iz=%d \n",ix,iy,iz);
	
	int indexcell=ix+(iy*maxcell)+(iz*maxcell*maxcell);

	//printf("ix=%d iy=%d iz=%d indexcell=%d\n",ix,iy,iz,indexcell);

	CELLLIST[indexcell][clist[indexcell]++]=ipart;

	if(clist[indexcell]>maxcellpart){

	  printf("clist[indexcell]=%d maxcellpart=%d\n",clist[indexcell],maxcellpart);	  printf("Something wrong: Increase maxcellpart size\n");
	  return 0;
	}

      }//ipart
      
      printf("indexing done\n");
      

       
       for(int ipart=0;ipart< Ncopy;ipart++){

	 if(ipart%smax==0)printf("fraction %d of 27 done\n",(ipart/smax));
	 
	 //for(int ineighbour=ipart+1;ineighbour< Ncopy;ineighbour++){

	int ix=img_xyz_0[ipart][0]/cellwidth;
	int iy=img_xyz_0[ipart][1]/cellwidth;
	int iz=img_xyz_0[ipart][2]/cellwidth;
	
	for(int celx=ix-1;celx<=ix+1;celx++){
	  int icelx=(celx+maxcell)%maxcell;
	  for(int cely=iy-1;cely<=iy+1;cely++){
	    int icely=(cely+maxcell)%maxcell;
	    for(int celz=iz-1;celz<=iz+1;celz++){
	      int icelz=(celz+maxcell)%maxcell;
	      
	      int indexcell=icelx+icely*maxcell+icelz*maxcell*maxcell;
	      
	      for(int icounter=0;icounter< clist[indexcell];icounter++){
		int ineighbour=CELLLIST[indexcell][icounter];
		if(ipart!=ineighbour){

		xr=img_xyz_0[ipart][0]-img_xyz_0[ineighbour][0];
		yr=img_xyz_0[ipart][1]-img_xyz_0[ineighbour][1];
		zr=img_xyz_0[ipart][2]-img_xyz_0[ineighbour][2];
		
		xr=xr-tildcopy*round(zr*iLcopy);
		xr=xr-Lcopy*round(xr*iLcopy);
		yr=yr-Lcopy*round(yr*iLcopy);
		zr=zr-Lcopy*round(zr*iLcopy);
		dr=sqrt(xr*xr+yr*yr+zr*zr);
		if(dr<1.4){
		  neighbourLIST[ipart][nlist[ipart]++]=ineighbour;
		  //neighbourLIST[ineighbour][nlist[ineighbour]++]=ipart;
		  
		}
		
		}//if ipart!=ineighbour	

	      }//icounter
	    }//cellz
	  }//celly
	} //cellx
	
	//}//ineighbour		

  	 
	
      }//ipart
      printf("Neighbourlist of image completed through cell list\n");
      /* //========================================================= */   
      




      


      
      ntot=0;ncl=0;
      for(i=0;i<Ncopy;i++)img_ID[i]=-1;
      for(i=0;i<=Ncopy;i++)img_ns[i]=0;
      //for(i=0;i<img_iplist;i++){
      for(i=0;i<Ncopy;i++){
	
      	iroot=i;
	
      	if(img_ID[iroot]==-1){
      	  ncl++;
      	  img_ID[iroot]=ncl;
      	  n=1;
      	  img_LIST[n]=iroot;
      	  k1=n;
      	  k2=k1;
      	  while(k1<=k2){
      	    for(k=k1;k<=k2;k++){
      	      iroot=img_LIST[k];

	      //for(j=0;j<img_iplist;j++){
	      for(j=0;j<nlist[iroot];j++){

		

		inew=neighbourLIST[iroot][j];
      		if(img_ID[inew]==-1){
		  
      		  //xr=img_xyz_0[iroot][0]-img_xyz_0[inew][0];
      		  //yr=img_xyz_0[iroot][1]-img_xyz_0[inew][1];
      		  //zr=img_xyz_0[iroot][2]-img_xyz_0[inew][2];
		  // 
      		  ////xr=xr-tildcopy*round(zr/Lcopy);
      		  //xr=xr-Lcopy*round(xr/Lcopy);
      		  //yr=yr-Lcopy*round(yr/Lcopy);
      		  //zr=zr-Lcopy*round(zr/Lcopy);
      		  //dr=sqrt(xr*xr+yr*yr+zr*zr);
		  
      		  //if(dr<1.4){
      		    img_ID[inew]=ncl;
      		    n++;
      		    img_LIST[n]=inew;
		    //}
      		}
      	      }//loop j
      	    }//loop k
      	    k1=k2+1;
      	    k2=n;
	    
	  }//while
      	  img_ns[n]++;
	  
      	  ntot+=n;
	  
	}//if
	
      	if(n>img_smax)img_smax=n;
	
      }//loop i
      

      printf("img_smax=%d smax=%d \n",img_smax,smax);
      Pdir=0;
      
      if(img_smax==3*smax){Pdir=1;printf("percolated one direction\n");}
      else if(img_smax==9*smax){Pdir=2;printf("percolated two directions\n");}
      else if(img_smax==27*smax){Pdir=3;printf("percolated three directions\n");
      }
      else {Pdir=0;printf("Not percolated\n");}
	
      
      
      

    }//if percolation
      
    
    //for(i=0;i<=N;i++)if(ns[i]>0)printf("%d %d\n",i,ns[i]);
      
    for(i=0;i<N;i++)if(ns[i]>0)tot_ns[i]+=ns[i];  
      
    sprintf(FILE,"outfiles/ncyc_smax_totactive_g0.069_s1.dat");
    fp=fopen(FILE,"a");
    fprintf(fp,"%d %d %d %d\n",cycle,smax,iplist,Pdir);
    fclose(fp);
    
    sprintf(FILE,"outfiles/Dist_ns_g0.069_s1_cycle%d.dat",cycle);
    fp=fopen(FILE,"w");
    for(i=0;i<=N;i++)if(tot_ns[i]>0)fprintf(fp,"%ld %d\n",i,tot_ns[i]);
    fclose(fp);

    /*
    sprintf(FILE,"outfiles/Smax_ns_g0.069_s1_cycle%d.xyz",cycle);
    fp=fopen(FILE,"w");
    fprintf(fp,"%d\n",smax);
    fprintf(fp,"Atoms\n");
    for(i=0;i<=smax;i++)fprintf(fp,"%f %f %f\n",xsmax[i][0],xsmax[i][1],xsmax[i][2]);
    fclose(fp);

    */
    
    //for(i=0;i<=N;i++)if(tot_ns[i]>0)printf("%d %d\n",i,tot_ns[i]);


//    int binsize=1000;
//    double binof=2.0;
//    double bincount[binsize],tot_count;
//    double b1,b2,bw;
//    for(i=0;i<binsize;i++)bincount[i]=0.0;
//    tot_count=0.0;
//    for(i=1;i<=N;i++){
//      index=log(i*1.0)/log(binof)+1;
//      bincount[index]+=tot_ns[i];
//      tot_count+=tot_ns[i];
//    }
//    
//    sprintf(FILE,"bin_Dist_ns_%3.2f.dat",gmax);
//    fp=fopen(FILE,"w");
//    for(i=0;i<binsize;i++){
//      b1=pow(binof,i-1);
//      b2=pow(binof,i)-1;
//      bw=sqrt(b1*b2);
//      if(bincount[i]>0)
//    	fprintf(fp,"%e %e\n",bw,bincount[i]/((b2-b1+1)*tot_count));
//    }
//    fclose(fp);
    
    
    
    free(xyz_0); 
    free(xyz_1); 
    free(PLIST); 
    free(LIST); 
    free(ID); 
    free(ns); 
    free(tot_ns); 
    
    return 0;
} //end main


