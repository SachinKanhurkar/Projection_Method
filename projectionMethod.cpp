#include <iostream>
#include <vector>
#include<stdio.h>
#include <stdlib.h>
#include <cstring>
#include <cmath>
#include <fstream>

using namespace std;
int i, j;

int main ()
{
int maxstep =200;
double mu= 0.1;
double rho =1;
int Iter= 10;
double beta =1;

int imax=16;int jmax=16;float Lx=1;float Ly=1;
 float dx=Lx/imax; float dy=Ly/jmax;
 cout<<"dx="<<dx<<endl;
 float  un=1;float us=0;float ve=0;float vw=0;
 float time_initial=0; float dt=0.005;
 int time_step;


float u[imax+1][jmax+2];float u_star[imax+1][jmax+2];
float v[imax+2][jmax+1];float v_star[imax+2][jmax+1];
float p[imax+2][jmax+2];float p_old[imax+2][jmax+2];
float c[imax+1][jmax+1];
float x[imax+1][jmax+1];float y[imax+1][jmax+1];

//initialize all array
for(i=0;i<=imax;i++)
for (j=0;j<=jmax+1;j++)
{u[i][j]=u_star[i][j]=0;}



for (i=0;i<=imax+1;i++)
for(j=0;j<=jmax;j++)
{v[i][j]=v_star[i][j]=0;}

for(i=0;i<=imax+1;i++)
for(j=0;j<=jmax+1;j++)
{p[i][j]=p_old[i][j]=0;
c[i][j]=1/(2/(dx*dx) + 2/(dy*dy));
}
   
	for (j=2;j<jmax;j++)
	{
	c[1][j]=1/(1/(dx*dx) + 2/(dy*dy));//left side boundary except corner
	}
	for(j=2;j<jmax;j++)
	{
        c[imax][j]=1/(1/(dx*dx) + 2/(dy*dy));//right side boundary except corner
	}
        for (i=2;i<imax;i++)
	{
	c[i][1]=1/(1/(dx*dx) + 2/(dy*dy));//south boundary except corner
	cout<<c[i][1]<<endl;
	}
	for (i=2;i<imax;i++) 
	{c[i][jmax]=1/(1/(dx*dx) + 2/(dy*dy));}//north boundary except corner
	
	c[1][1]=1/(1/(dx*dx) + 1/(dy*dy)); //4 corner point 
	c[1][jmax]=1/(1/(dx*dx) + 1/(dy*dy));
	c[imax][1]=1/(1/(dx*dx) + 1/(dy*dy));
	c[imax][jmax]=1/(1/(dx*dx) + 1/(dy*dy));

	for (i=0;i<=imax;i++)
	{for (j=0;j<=jmax;j++)
        {
	x[i][j]=dx*(i);y[i][j]=dy*(j);
	cout<<x[i][j]<< "  ";
	}
	cout<< endl;
	}
	
/************************Start time loop*****************************************************************/
for (int t=1;t<=maxstep;t++)
 { //boundary coondition for velocity ; We need the tengential velocoty at the each boundary but not the normal velocity ,which is given .
      	      	cout<<"Time "<< t << endl;
      for (i=0;i<=imax;i++)
      	{
      	u[i][0]=2*us-u[i][1];
      	
      	}
          
	 for (i=0;i<=imax;i++)
      	{u[i][jmax+1]=2*un-u[i][jmax];}
    
	 for (j=0;j<=jmax;j++)
     	 {v[0][j]=2*vw-v[1][j];}
	
	for (j=0;j<=jmax;j++)
      	{v[imax+1][j]=2*ve-v[imax][j];}
// temporary u-velocity
        for(i=1;i<imax;i++)
	for(j=1;j<=jmax;j++)    
       {       
       u_star[i][j]=u[i][j]+dt*((-0.25/(dx))*((u[i+1][j]+u[i][j])*(u[i+1][j]+u[i][j])-(u[i][j]+u[i-1][j])*(u[i][j]+u[i-1][j])+(u[i][j+1]+u[i][j])*(v[i+1][j]+v[i][j])-(u[i][j]+u[i][j-1])*(v[i+1][j-1]+v[i][j-1]))+(mu/(dx*dx))*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]-4*u[i][j]));
	}
// temporary v-velocity	
    for( i=1;i<=imax;i++)
    for( j=1;j<jmax;j++ )      
     {
      v_star[i][j]=v[i][j]+dt*((-0.25/dx)*((u[i][j+1]+u[i][j])*(v[i+1][j]+v[i][j])-(u[i-1][j+1]+u[i-1][j])*(v[i][j]+v[i-1][j])+(v[i][j+1]+v[i][j])*(v[i][j+1]+v[i][j])-(v[i][j]+v[i][j-1])*(v[i][j]+v[i][j-1]))+(mu/(dx*dx))*(v[i+1][j]+v[i-1][j]+v[i][j+1]+v[i][j-1]-4*v[i][j]));
     }
cout<<"test5"<<endl;
//*****************************solve the pressure poission equation*********************************************************************//
	for (int k=1;k<=100;k++)
	{ cout<<"iteration:"<<k<<endl;
	for (i=1;i<=imax;i++)
	for (j=1;j<=jmax;j++)
	{
	p[i][j]=beta*c[i][j]*((p[i+1][j]+p[i-1][j]+p[i][j+1]+p[i][j-1])/(dx*dx)-(rho/(dx*dt))*(u_star[i][j]-u_star[i-1][j]+v_star[i][j]-v_star[i][j-1]))+(1-beta)*p[i][j];
	
	}
	
	}
        for (i=0;i<=imax+1;i++)
	{for (j=0;j<=jmax+1;j++)
	{
	cout<<p[i][j]<< "/";
	
	}
	}
//************************************************correct the velocity******************************************************************//
for (i=1;i<imax;i++)
for (j=1;j<=jmax;j++)
{
u[i][j]=u_star[i][j]-(dt/(dx*rho))*(p[i+1][j]-p[i][j]);
cout<<u[i][j]<<endl;
}

for (i=1;i<=imax;i++)
for (j=1;j<jmax;j++)
{
v[i][j]=v_star[i][j]-(dt/(dy*rho))*(p[i][j+1]-p[i][j]);
}
cout<<"test4"<<endl;
}

for (i=0;i<=imax;i++)
{
for (j=0;j<=jmax+1;j++)
{

cout<<u[i][j]<<"/";
}
cout<<"newrow_u"<<endl;
}
cout<<endl;

for (i=0;i<=imax+1;i++)
{
for (j=0;j<=jmax;j++)
{

cout<<v[i][j]<<"/";
}
cout<<"new row_v"<<endl;
}
 /***************************plotting result*************************************************/
 float uu[imax+1][jmax+1];float vv[imax+1][jmax+1];float w[imax+1][jmax+1];float u_v[imax+1][jmax+1];
 for (i=0;i<=imax;i++)
 {for (j=0;j<=jmax;j++)
 {
 uu[i][j]=0.5*(u[i][j+1]+u[i][j]);
 vv[i][j]=0.5*(v[i+1][j]+v[i][j]);
 w[i][j]=(u[i][j+1]-u[i][j]-v[i+1][j]+v[i][j])/(2*dx);
u_v[i][j]=sqrt(uu[i][j]*uu[i][j]+vv[i][j]*vv[i][j]);
 }
 }
  
 ofstream out1;
 out1.open("vorticity.txt"); 
 for (j=0;j<=jmax;j++)
 {for (i=0;i<=imax;i++)
 {
 out1<<w[i][j]<<endl;
  }
  out1<<endl;
  }
  
  ofstream out2;
 out2.open("u_velocty.txt"); 
 for (j=0;j<=jmax;j++)
 {for (i=0;i<=imax;i++)
 {
 out2<<uu[i][j]<<endl;
  }
  out2<<endl;
  }
  
  ofstream out3;
 out3.open("v_velocty.txt"); 
 for (j=0;j<=jmax;j++)
 {for (i=0;i<=imax;i++)
 {
 out3<<vv[i][j]<<endl;
  }
  out3<<endl;
  }

  ofstream out4;
 out4.open("Resultant_velocty.txt"); 
 for (j=0;j<=jmax;j++)
 {for (i=0;i<=imax;i++)
 {
  out4<<u_v[i][j]<<endl;
  }
  out4<<endl;
  }
return 0;
}

