#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define H 1
#define W 1
#define nx 129
#define ny 129
#define Re 400
main()
{
	int i,j,tmstep=0;
	double dx,dy,dt,beta,errp,erru,sumu,sump;
	double u[nx+2][ny+2],u_new[nx+2][ny+2],u2[nx+2][ny+2],uv[nx+2][ny+2];
	double v[nx+2][ny+2],v2[nx+2][ny+2],sf[nx+2][ny+2];
	double p[nx+2][ny+2],p_new[nx+2][ny+2];
	double RHSU[nx+2][ny+2],RHSV[nx+2][ny+2];
		
	dx = (double)W/(nx-1); 
	dy = (double)H/(ny-1);
	dt = 0.001;
	beta = dx/dy;
	
	//INITIAL
	for(i=1;i<=nx+1;i++){
		for(j=1;j<=ny+1;j++){
			p[i][j] = 0;
			p_new[i][j] = 0;
			
			u[i][j] = 0;
			u_new[i][j] = 0;
			RHSU[i][j] = 0;
			v[i][j] = 0;
			RHSV[i][j] = 0;
			uv[i][j] = 0;
			sf[i][j] = 0;
		
		}			
	}

	erru=1;
	while(erru>0.00001){
		
		for(i=2;i<=nx-1;i++){
			u[i][1] = -u[i][2];
			u[i][ny+1] = 2-u[i][ny];
		}		
		for(j=2;j<=ny-1;j++){
			v[1][j] = -v[2][j];
			v[nx+1][j] = -v[nx][j];
		}
		
		for(i=2;i<=nx;i++){
			for(j=2;j<=ny;j++){
				u2[i][j] = 0.25*(u[i][j]+u[i-1][j])*(u[i][j]+u[i-1][j]); 	
				v2[i][j] = 0.25*(v[i][j]+v[i][j-1])*(v[i][j]+v[i][j-1]); 
	
			}
		}
	
		for(i=2;i<=nx-1;i++){
			for(j=2;j<=ny-1;j++){
				uv[i][j] = 0.25*(u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j]);
			}
		}
	
		for(i=2;i<=nx-1;i++){
			for(j=2;j<=ny;j++){
				RHSU[i][j] = u[i][j] - (dt/dx)*(u2[i+1][j]-u2[i][j]) - (dt/dy)*(uv[i][j]-uv[i][j-1]) 
				+ dt*(u[i+1][j]-2*u[i][j]+u[i-1][j])/(dx*dx*Re)+(dt/(Re*dy*dy))*(u[i][j-1]-2*u[i][j]+u[i][j+1]);	
			}
		}
	
		for(i=2;i<=nx;i++){
			for(j=2;j<=ny-1;j++){
				RHSV[i][j] = v[i][j] - ((double)dt/dy)*(v2[i][j+1]-v2[i][j]) - ((double)dt/dx)*(uv[i][j]-uv[i-1][j])+((double)dt/(dx*dx*Re))*(v[i+1][j]-2*v[i][j]+v[i-1][j])+((double)dt/(Re*dy*dy))*(v[i][j-1]-2*v[i][j]+v[i][j+1]);	
			}
		}
		
		
		errp = 1;
		while(errp>0.00001){
			for(i=1;i<=nx+1;i++){
				for(j=1;j<=ny+1;j++){	
					p[1][j] = p[2][j];
					p[nx+1][j] = p[nx][j];
					p[i][1] = p[i][2];
					p[i][ny+1] = p[i][ny];					
				}
			}
			errp = 0;
			sump = 0;
			for(i=2;i<=nx;i++){
				for(j=2;j<=ny;j++){
					p_new[i][j] = (double)(p[i+1][j]+p[i-1][j] + beta*beta*(p[i][j+1]+p[i][j-1])) - (double)(dx*dx/dt)*((RHSU[i][j]-RHSU[i-1][j])/dx + (RHSV[i][j]-RHSV[i][j-1])/dy);
					p_new[i][j] = 0.5*p_new[i][j]/(1+beta*beta);
					sump += fabs(p_new[i][j]);
					errp += fabs(p[i][j]-p_new[i][j]);
					p[i][j] = p_new[i][j];							
				}
			}
			
			for(i=1;i<=nx+1;i++){
				for(j=1;j<=ny+1;j++){	
					p[1][j] = p[2][j];
					p[nx+1][j] = p[nx][j];
					p[i][1] = p[i][2];
					p[i][ny+1] = p[i][ny];					
				}
			}
			errp = (double)errp/sump;
		//	printf("%lf	",errp);		
		}
		//printf("\n");	
		
		//UPDATING VELOCITIES
		sumu = 0;
		erru = 0;
		for(i=2;i<=nx-1;i++){
			for(j=2;j<=ny;j++){
				u_new[i][j] = -(dt/dx)*(p[i+1][j]-p[i][j]) + RHSU[i][j];
				sumu += fabs(u_new[i][j]);
				erru += fabs(u_new[i][j]-u[i][j]);
				u[i][j]=u_new[i][j];
			}
		}
		for(i=2;i<=nx;i++){
			for(j=2;j<=ny-1;j++){
				v[i][j] = -(dt/dy)*(p[i][j+1]-p[i][j]) + RHSV[i][j];
			}
		}			
		erru = (double)erru/sumu;	
	//	printf("%lf\n",erru);
	}
		
		
	//RESULTS PRINTING
		FILE *fa,*fb,*fc,*fd,*fe,*ff;
	
	fa = fopen("MACpressure.txt", "w+");   
		fprintf(fa,"pressure:\n\n");
		for(j=ny;j>=2;j--){
			for(i=2;i<=nx;i++){
				fprintf(fa,"%lf	",p[i][j]);
			}
			fprintf(fa,"\n");
		}
		fclose(fa);
	
	fb = fopen("MACu_velocities.txt", "w+");  	
		fprintf(fb,"u velocities:\n\n");
		for(j=ny;j>=1;j--){
			for(i=1;i<=nx;i++){
				fprintf(fb,"%lf	",u[i][j]);
			}
			fprintf(fb,"\n");
		}
	
	fc = fopen("MACv_velocities.txt", "w+");  	
		fprintf(fc,"v velocities:\n\n");	
		for(j=ny;j>=1;j--){
			for(i=1;i<=nx;i++){
				fprintf(fc,"%lf	",v[i][j]);
			}
			fprintf(fc,"\n");
		}
		
	fd = fopen("MAC_stream.txt", "w+");
		fprintf(fd,"Stream function values:\n\n");
		for(i=2;i<=nx-1;i++){
			for(j=2;j<=ny;j++){
				sf[i][j] = sf[i-1][j] + v[i][j]*dx;
				fprintf(fd,"%lf	",sf[i][j]);
			}
			fprintf(fd,"\n");
		}
		
	fe = fopen("mid u velocities.txt", "w+");  	
		fprintf(fe,"mid u velocities:\n\n");	
		for(j=ny;j>=1;j--){
				fprintf(fe,"%lf	\n",u[65][j]);
		}	
		
	ff = fopen("mid v velocities.txt", "w+");  	
		fprintf(ff,"mid v velocities:\n\n");	
			for(i=1;i<=nx;i++){
				fprintf(ff,"%lf	\n",v[i][65]);
			}
	printf("\n\nResults have been printed in the output text files\n\n");		
}
	



