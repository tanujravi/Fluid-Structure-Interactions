#include<iostream>
#include<fstream>
#include<cmath>
#include <sys/types.h>
#include <sys/stat.h>

using namespace std;

const double Re=10.0; //Reynolds Number
const double length1= 6.0;//Length
const double length2= 1.0;//Height
const double length3= 0.4;//Height
const double w = 4.0;
const double Lc = 2.0;
const double Ls = 1.0;
const double co = 0.125/2;
double dt= 0.0001;// TimeStep
const int N= 151;
const int M= 51;
const int NLagMem = 350;
const int NLagIn = 75;
const int Iin = (1 - w/length1)*(N-1)/2;
//const int Iout = (1 + w/length1)*(N-1)/2;
const int Iout = 125;
const int Jin = (1 - length3/length2)*(M-1)/2 ;
const int Jout = (1 + length3/length2)*(M-1)/2 - 1;
const double rho= 1.0;
double tolerance=pow(10,-6);
const double NT=3000000;//Max number of iterations allowed
const double dsMem = w/(NLagMem - 1);
const double dsIn = length3/(NLagIn - 1);
//feedback controls
double alpha = 100;
double beta = 1.2;
double P_blood = 1;
double Ui_centre = 1.0;
double Uo_centre;
//const double Ks = 0*pow(10,12);
const double Ks =75000;
const double Kb = 0.0000;
// defining variables
double u_infinity= 1.0;

double fx[N+1][M+1], fx_n[N+1][M+1], fx_h[N+1][M+1];
double fy[N+1][M+1], fy_n[N+1][M+1], fy_h[N+1][M+1];
double p[N+1][M+1],p_error[N+1][M+1]={0.0};
double u[N+1][M+1], un[N+1][M+1],u_n_1[N+1][M+1]={0.0};
double v[N+1][M+1], vn[N][M+1],v_n_1[N+1][M+1]={0.0};
double l[N+1],b[N+1],c[N+1],d[N+1],x[N+1], xm[M+1];//TDMA Coefficients

//Momentum Equation Coefficient Functions
double a(int ix, int iy,int z);
double a_w(int ix, int iy,int z);
double a_e(int ix, int iy,int z);
double a_s(int ix, int iy,int z);
double a_n(int ix, int iy,int z);

//Pressure Correction Equation Coefficients
double p_e(int ix, int iy);
double p_s(int ix, int iy);
double p_n(int ix, int iy);
double p_w(int ix, int iy);
double dx= length1/(N-1);
double dy= length2/(M-1);

//double viscosity= (rho*u_infinity*length2)/rey ;
//double D1= viscosity/dx;
//double D2= viscosity/dy;

double R1=dy/dx;
double R2=dx/dy;

//Momentum-Equation Solver
void gs_u(void);
void gs_v(void);

//Pressure Correction Solver
void gs_p(void);

void TDMAm(void);//TDMA X-Sweep_Momentum Eqns
void TDMAn(void);//TDMA Y-Sweep_Momentum Eqns
void TDMApm(void);
void SIMPLE(int mn, double t);

class Lagrange{
public:
	double xLag, yLag, xs, ys, xh, yh, Zx, Zy;
	double ULag, VLag;
	double Fx_Lag, Fx_Lag_next, Fhx_Lag;
	double Fy_Lag, Fy_Lag_next, Fhy_Lag;
	int positionx, positiony;


	//float ds;


};

double d_s(Lagrange A, Lagrange B)
    {
        double ds;
        ds = sqrt((A.xLag - B.xLag)*(A.xLag - B.xLag) + (A.yLag - B.yLag)*(A.yLag - B.yLag));
        return ds;

    }
double phi(double r1)
{

	double r;
    if(r1 < -1 && r1 >= -2)
        r = (5 + 2*r1 - sqrt(-7 - 12*r1- 4*r1*r1))/8;
    else if(r1 >= -1 && r1 < 0)
        r = (3 + 2*r1 + sqrt(1 - 4*r1- 4*r1*r1))/8;
    else if(r1 >= 0 && r1 < 1)
        r = (3 - 2*r1 + sqrt(1 + 4*r1- 4*r1*r1))/8;
    else if(r1 >= 1 && r1 < 2)
        r = (5 - 2*r1 - sqrt(-7 + 12*r1- 4*r1*r1))/8;
    else
        r = 0;
	/*if(abs(r1) <= 2)
		r = (1 + cos(3.14*r1/2))/4;
	else
		r = 0;*/
    //cout<<"r = "<<r<<"  "<<"r1 = "<<r1<<endl;
	return r;
}

double deltah(double x, double y, int k)
{
	double d_h;

	d_h = phi(x/dx)*phi(y/dy)/(dx*dy);
    //cout<<"d_h = "<<d_h<<"   "<<x<<endl;
	return d_h;

}

int del_kro(int m, int l)
{
	if(m == l)
		return 1;
	else
		return 0;
}

double UpdateMembranePosition1(Lagrange L, int k, int rq);
double UpdateMembranePosition(Lagrange L, int k, int rq);
double ForceCalculate(Lagrange L, Lagrange Lp, Lagrange Ln, Lagrange**A, int k);
double ForceCalculate_h(Lagrange L, Lagrange Lp, Lagrange Ln, Lagrange**A, int k);
double EulerianForce(Lagrange**A, int p, int q, int k);
double EulerianForce_h(Lagrange*A, int p, int q, int k);
double LagrangianVelocity(Lagrange L, int k);

int main()
	{

		/*if (mkdir("./TimeData", ACCESSPERMS)!=0);
			cout<<"Couldn't Create Directory"<<endl;*/

		ofstream fout30, fout31, fout32, fout33;
        fout30.open("Lagrangian Pts.dat");
        fout31.open("Lagrangian Pts1.dat");
        fout32.open("EulerainForce1.dat");
        fout33.open("Data.dat");
        Lagrange** bdry = new Lagrange*[4];
        for(int i = 0; i < 4; i++)
        	bdry[i] = new Lagrange[NLagMem];



              int i=0;
        int j=0;

        for(int i = 0; i < 2; i++)
        {
        	for(int j = 0; j < NLagMem; j++)
        	{
        		bdry[i][j].positiony = i;
        		bdry[i][j].positionx = j;
        		//bdry[i][j].xLag = j * dsMem + (length1-w)/2;
        		//bdry[i][j].yLag = i * length3 + (length2-length3)/2;

                bdry[i][j].xLag = j * dsMem + (length1-w)/2;
        		if(abs(j * dsMem-Lc)>=Ls/2)
        		bdry[i][j].yLag = i * length3 + (length2-length3)/2;
                if(abs(j * dsMem-Lc)<Ls/2)
                {
                if(i==0)
                bdry[i][j].yLag = (length2-length3)/2 + (co/2)*(1+cos(2*2*3.14*(j * dsMem-Lc)*(1.0/Lc))) ;
                if(i==1)
                bdry[i][j].yLag = (length2+length3)/2 - (co/2)*(1+cos(2*2*3.14*(j * dsMem-Lc)/Lc)) ;}

        		bdry[i][j].Zx = j * dsMem + (length1-w)/2;
        		if(abs(j * dsMem-Lc)>=Ls/2)
        		bdry[i][j].Zy = i * length3 + (length2-length3)/2;
                if(abs(j * dsMem-Lc)<Ls/2)
                {
                if(i==0)
                bdry[i][j].Zy = (length2-length3)/2 + (co/2)*(1+cos(2*2*3.14*(j * dsMem-Lc)*(1.0/Lc))) ;
                if(i==1)
                bdry[i][j].Zy = (length2+length3)/2 - (co/2)*(1+cos(2*2*3.14*(j * dsMem-Lc)/Lc)) ;}

        		bdry[i][j].xh = j * dsMem;
        		bdry[i][j].yh = i*length2*0.5;
        		bdry[i][j].Fx_Lag = 0;
        		bdry[i][j].Fhx_Lag = 0;
        		bdry[i][j].Fy_Lag = 0;
        		bdry[i][j].Fhy_Lag = 0;
        		bdry[i][j].ULag = 0;
        		bdry[i][j].VLag = 0;

        		fout30<<bdry[i][j].Zx<<" "<<bdry[i][j].Zy<<endl;
        	}
        }

        for(int i = 2; i < 4; i++)
        {
        	for(int j = 0; j < NLagIn; j++)
        	{
        		bdry[i][j].positiony = i;
        		bdry[i][j].positionx = j;
        		if(i==2)
                bdry[i][j].xLag = (length1-w)/2;
                if(i==3)
                bdry[i][j].xLag = w + (length1-w)/2;
        		bdry[i][j].yLag = j * dsIn + (length2-length3)/2;
        		bdry[i][j].xh = j * dsMem;
        		bdry[i][j].yh = i*length2*0.5;
        		bdry[i][j].Fx_Lag = 0;
        		bdry[i][j].Fhx_Lag = 0;
        		bdry[i][j].Fy_Lag = 0;
        		bdry[i][j].Fhy_Lag = 0;
        		bdry[i][j].ULag = 0;
        		bdry[i][j].VLag = 0;

        		//fout30<<bdry[i][j].xLag<<" "<<bdry[i][j].yLag<<endl;
        	}
        }

        //cout<<3;

        cout<<"Iin = "<<Iin<<"Iout = "<<Iout<<"Jin = "<<Jin<<"Jout = "<<Jout<<endl;
        double t = 0;
        double z = 0;
        double maxu = 0, erru = 0;
        double max = 1;
        ofstream fout20;
        fout20.open("PressDiff.dat");
        //Boundary Conditions-Velocity
                for (int i=M*(length2 - length3)*0.5;i<M*(length2 + length3)*0.5;i++)
                {
                    u[0][i]=0;
                }

                for (int i=1;i<M;i++)
                {
                    //u[0][i]=1;//Left_Velocity Inlet
                    u[N-1][i]=0;//Right
                }
                for (int i=1;i<N;i++)
                {
                    v[i][0]=0.0;//Bottom
                    v[i][M-1]=0.0;//Top
                }
                for (int i=0;i<=M-1;i++)
                {
                    v[0][i]=-v[1][i];//Left
                    v[N][i]=-v[N-1][i];//Right
                }
                for (int i=0;i<=N-1;i++)
                {
                    u[i][0]=-u[i][1];//Bottom
                    u[i][M]=-u[i][M-1];//Top
                }
                /*p[0][0] = 0;
                p[N][0] = 0;
                p[0][M] = 0;
                p[N][M] = 0;*/

        while(max > 1e-6)//1e-6 is convergence criteria for steady state

        {
            maxu = 0;
            z++;
            t =z*dt;
            cout<<"Time = "<<t<<"s"<<endl<<endl;
            char time_data[50], Bdry_data[50], time_datav[50], eulerian_force[50];

            // Updating previous timestep values-u velocity
        	for(int ix=0;ix<N;ix++)
        	{
            	for(int iy=0;iy<=M;iy++)
            	{
                	u_n_1[ix][iy]=u[ix][iy];//u_n_1 is velocity of previous timestep
            	}
       		}

       		// Updating previous timestep values-v velocity
       		for(int ix=0;ix<=N;ix++)
       		{
        		for(int iy=0;iy<M;iy++)
        		{
            		v_n_1[ix][iy]=v[ix][iy];
        		}
       		}

            for(int i = 0; i < N; i++)
            {
                for(int j = 0; j < M; j++)
                {
                    fx_n[i][j] = fx_h[i][j];
                    fy_n[i][j] = fy_h[i][j];
                }
            }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 //                       SIMPLE algorithm
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	 	/*for(int i = 0; i < 2; i++)
        {
        	for(int j = 1; j < NLag - 1; j++)
        	{*/
            for(int i = 0; i < 2; i++)
        {
        	for(int j = 0; j <= NLagMem -1; j++)
        	{
                //cout<<bdry[i][j].xLag<<endl;
        		//bdry[i][j].xLag = UpdateMembranePosition(bdry[i][j], 1, 2);
                bdry[i][j].yLag = UpdateMembranePosition(bdry[i][j], 2, 2);
                //bdry[i][j].yLag = PeriodicUpdate(bdry[i][j], 2);*/
        		//fout30<<bdry[i][j].xLag<<" "<<bdry[i][j].yLag<<endl;


        	}
        }

         for(int i = 0; i < 2; i++)
        {
        	for(int j = 1; j < NLagMem-1; j++)
        	{
        	    bdry[i][j].Fx_Lag = ForceCalculate(bdry[i][j], bdry[i][j+1], bdry[i][j-1], bdry, 1);
        		bdry[i][j].Fy_Lag = ForceCalculate(bdry[i][j], bdry[i][j+1], bdry[i][j-1], bdry, 2);
        		//fout31<<bdry[i][j].xLag<<" "<<bdry[i][j].yLag<<"  "<<bdry[i][j].Fx_Lag<<"  "<<bdry[i][j].Fy_Lag<<endl;

        	}
        }

    double P_in = 0;
    double P_ex = 0;
    double Q_in = 0;
    double Q_out;
    double sumg = 0;
    double g[Jout-Jin+1];
    int b;
    b = length3*(M-1) + 1;
    double R[b], I[b];
    int r = 0;
    int im = 0;
    ifstream fin1, fin2;
    ofstream fout1;

    fin1.open("Real Part.txt");
    fin2.open("Imaginary Part.txt");

    while(!fin1.eof())
    {
        fin1>>R[r];
        r++;
    }

    while(!fin2.eof())
    {
        fin2>>I[im];
        im++;
    }
    cout<<"---------------------------------"<<endl;

    for(int j = 0; j <= Jout-Jin; j++)
         cout<<R[j]<<endl;

cout<<"---------------------------------"<<endl;
    for(int j = 0; j <= Jout-Jin; j++)
         cout<<I[j]<<endl;
    cout<<"---------------------------------"<<endl;

    r = 0;
    im = 0;
    fout1.open("Vel Profile.dat");
    for(int j = 0; j <= Jout-Jin; j++)
    {

        g[j] = -(R[r]*sin(2*3.142857*t) + I[im]*cos(2*3.142857*t));
        r++;
        im++;
        fout1<<u[0][i]<<endl;
    }


    for(int i = 0; i <= N; i++)
            {
        for(int j = 0; j <= M; j++)
         {
             if(j*dy > (length2 - length3)/2 && j*dy < (length2 + length3)/2 && i*dx > (length1-w)/2 && i*dx < (length1+w)/2)
               {P_in = P_in + p[i][j]*dy/(length3*(Iout-Iin));}
             else if (i*dx <= (length1-w)/2 || i*dx >= (length1+w)/2)
               {P_ex = P_ex + p[i][j]*dy/(2*Iin*length3);}
             else
               {P_ex = P_ex + p[i][j]*dy/((length2-length3)*(Iout-Iin));}
         }}


     /*for(int j = 0; j <= Jout-Jin; j++)
        g[j] = (4*(0.5*dy+j*dy)/length3 - 4*(0.5*dy+j*dy)*(0.5*dy+j*dy)/(length3*length3));*/

     for(int j = 0; j <= Jout-Jin; j++)
        Q_in = Q_in + g[j]*dy*Ui_centre;

     for(int j = 0; j <= Jout-Jin; j++)
        sumg = sumg + g[j]*dy;



        Q_out = Q_in - beta*(P_blood - (P_in - P_ex));

       cout<<"Pin = "<<P_in<<endl;
       cout<<"Pex = "<<P_ex<<endl;
       cout<<"Pin - Pex = "<<P_in - P_ex<<endl;
       Uo_centre = Q_out/sumg;
          //Uo_centre = Ui_centre;



        //Eulerian Force

        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < M; j++)
            {
                fx[i][j] = EulerianForce(bdry, i, j, 1);
                fy[i][j] = EulerianForce(bdry, i, j, 2);
                //fout32<<fx[i][j]<<"  "<<fy[i][j]<<endl;
            }
        }

       int cou = 0;
        for(int j = Jin; j <= Jout; j++)
            {
        fx[Iin][j] = alpha*(Ui_centre*g[cou] - u[Iin][j]);
        fy[Iin][j] = alpha*(0 - v[Iin][j]);
        cou++;}

        cout<<"F Inlet"<<endl;
        for(int j = Jin; j <= Jout; j++)
            {cout<<fx[Iin][j]<<endl; }


        cout<<"Uo_centre = "<<Uo_centre<<endl;
        cout<<"U_outlet"<<endl;

        cou = 0;
        for(int j = Jin; j <= Jout; j++)
            {
        cout<<u[Iout][j]<<endl;
        }
        cout<<"-------Diff Inlet------------"<<endl;
        cou = 0;
        for(int j = Jin; j <= Jout; j++)
            {
        cout<<Ui_centre*g[cou] - u[Iin][j]<<endl;
        cou++;
        }


        cou=0;
        for(int j = Jin; j <= Jout; j++)
            {
        fx[Iout][j] = alpha*(Uo_centre*g[cou] - u[Iout][j]);
        fy[Iout][j] = alpha*(0 - v[Iout][j]);
        cou++;}

        cout<<"-------Diff Outlet------------"<<endl;
        cou = 0;
        for(int j = Jin; j <= Jout; j++)
            {
        cout<<Uo_centre*g[cou] - u[Iout][j]<<endl;
        cou++;
        }

         cout<<"F outlet"<<endl;
        for(int j = Jin; j <= Jout; j++)
            {cout<<fx[Iout][j]<<endl; }



         for(int i = 0; i < N; i++)
            {
                for(int j = 0; j < M; j++)
                {
                    fx_n[i][j] = fx[i][j];
                    fy_n[i][j] = fy[i][j];
                }
            }

       SIMPLE(2, t);

       //double Nf = t/0.001;
           //if(floor(Nf)==ceil(Nf))
          //  {
                cout<<"Writing values to file."<<endl;
                fout20<<t<<" "<<p[Iout][M/2]-p[Iin][M/2]<<endl;
          //write values to a file after one time step converges
            int kt = round(t*1000);
            snprintf(time_data,sizeof time_data, "./TimeData/Data_%.2f.txt",t);
            snprintf(time_datav,sizeof time_datav, "./TimeDataV/%d.txt",kt);
            snprintf(eulerian_force,sizeof eulerian_force, "./eulerian_force/Data_%.2f.txt",t);
            ofstream ofp1(time_data, ios::out);
            ofstream ofp3(time_datav, ios::out);
            ofstream ofp4(eulerian_force, ios::out);
            int m;
           	for(int i=0;i<M;i++)
           	{
                m=25;
                ofp1<<0.5*(u[10][i]+u[10][i+1])/u_infinity<<" ";
                while (m<=N)
                {
                    ofp1<<(u[m][i]+u[m][i+1])*0.5/u_infinity<<" ";
                    m=m+25;
                }
                ofp1<<(dy*i)/length2<<endl;
            }

           for(int i = 0; i < N; i++){
			for(int j = 0; j < M; j++)
			{
				ofp3<<i*dx<<" "<<j*dy<<" "<<0.5*(u[i][j]+u[i][j+1])<<" "<<0.5*(v[i][j]+v[i+1][j])<<" "<<0.25*(p[i][j]+p[i+1][j]+p[i+1][j+1]+p[i][j+1])<<endl;
			}
		}

			//fout15.close();
            snprintf(Bdry_data,sizeof Bdry_data, "./BdryData/Data_%.3f.txt",t);
            ofstream ofp2(Bdry_data, ios::out);
            for(int i = 0; i < 2; i++)
            {
                for(int j = 0; j < NLagMem; j++)
                {
                    ofp2<<bdry[i][j].xLag<<" "<<bdry[i][j].yLag<<"  "<<bdry[i][j].Fhx_Lag<<"  "<<bdry[i][j].Fhy_Lag<<endl;
                }
            }
            for(int ix=0;ix<=N;ix++)
			{
           		for(int iy=0;iy<=M;iy++)
           		{
            		ofp4<<(ix)*dx+0.5*dx<<" "<<(iy)*dy+0.5*dy<<" "<<fx[ix][iy]<<" "<<fy[ix][iy]<<endl;
           		}
       		}

                //fout15.close();
			for(int ix=0;ix<N;ix++)
			{
           		for(int iy=0;iy<M;iy++)
           		{
            		erru = abs(u[ix][iy] - u_n_1[ix][iy]);
            		if(erru > maxu)
                	maxu =erru;
           		}
       		}
       		max = maxu;
    		cout<<"Velocity Error at the end of timestep: "<<max<<endl;
    		ofstream fout25;
    		fout25.open("Velocity Contour.dat");
    		//fout20<<t<<"  "<<u[N/2][M/2]<<endl;
    		for(int ix=1;ix<N+1;ix++)
			{
           		for(int iy=1;iy<M+1;iy++)
           		{
            		fout25<<(ix-1)*dx<<" "<<(iy-1)*dy<<" "<<0.5*(u[ix-1][iy]+u[ix-1][iy-1])<<" "<<0.5*(v[ix-1][iy-1]+v[ix][iy-1])<<endl;
           		}
       		}
       		//}
    	}
	}


// momentum equation coefficients

double a_n(int ix, int iy,int z){

    double ans=0.0;

    if(z==1)
     ans= dx*0.25*(v_n_1[ix][iy]+v_n_1[ix+1][iy])-R2/Re;
     else if(z==2)
     ans= dx*0.25*(v_n_1[ix][iy]+v_n_1[ix][iy+1])-R2/Re;

    return ans;
}

double a_s(int ix, int iy,int z){

    double ans=0.0;

        if(z==1)
        ans= -dx*0.25*(v_n_1[ix][iy-1]+v_n_1[ix+1][iy-1])-R2/Re;
        else if(z==2)
        ans= -dx*0.25*(v_n_1[ix][iy-1]+v_n_1[ix][iy+1])-R2/Re;

    return ans;
}


double a_e(int ix, int iy,int z){

    double ans=0.0;

    if(z==1)
      ans= dy*0.25*(u_n_1[ix][iy]+u_n_1[ix+1][iy])-R1/Re;
      else if(z==2)
      ans= dy*0.25*(u_n_1[ix][iy]+u_n_1[ix][iy+1])-R1/Re;

    return ans;
}

double a_w(int ix, int iy,int z){

    double ans=0.0;

    if(z==1)
    ans= -dy*0.25*(u_n_1[ix][iy]+u_n_1[ix-1][iy])-R1/Re;
    else if(z==2)
    ans= -dy*0.25*(u_n_1[ix-1][iy]+u_n_1[ix-1][iy+1])-R1/Re;


    return ans;
}


double a(int ix, int iy,int z){

    double ans=0.0;

    if(z==1)
     //ans= dx*dy/dt +(a_w(ix,iy,1)+a_e(ix,iy,1)+a_n(ix,iy,1)+a_s(ix,iy,1))+4.0*(R1/Re+R2/Re);
     ans= (dx*dy/dt)+(a_w(ix,iy,1)+a_e(ix,iy,1)+a_n(ix,iy,1)+a_s(ix,iy,1))+4.0*(R1/Re+R2/Re);
      //ans= (a_w(ix,iy,1)+a_e(ix,iy,1)+a_n(ix,iy,1)+a_s(ix,iy,1));
     else if(z==2)
     //ans= dx*dy/dt +(a_w(ix,iy,2)+a_e(ix,iy,2)+a_n(ix,iy,2)+a_s(ix,iy,2))+4.0*(R1/Re+R2/Re);
    ans= (dx*dy/dt)+(a_w(ix,iy,2)+a_e(ix,iy,2)+a_n(ix,iy,2)+a_s(ix,iy,2))+4.0*(R1/Re+R2/Re);
    //ans= (a_w(ix,iy,2)+a_e(ix,iy,2)+a_n(ix,iy,2)+a_s(ix,iy,2));
    return ans;
}

// pressure correction equation coefficients

double p_e(int ix, int iy){

    double answer=0.0;

    answer= -1.0*dy*dy*1.0/a(ix,iy,1);

    return answer;

}

double p_n(int ix, int iy){

    double answer=0.0;


    answer= -1.0*dx*dx*1.0/a(ix,iy,2);


    return answer;

}

double p_w(int ix, int iy){

    double answer=0.0;

    answer= -1.0*dy*dy*1.0/a(ix-1,iy,1);

    return answer;

}

double p_s(int ix, int iy){

    double answer=0.0;

    answer= -1.0*dx*dx*1.0/a(ix,iy-1,2);


    return answer;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void gs_u(void){
//Y-Sweep Substitution
for(int iy=1;iy<M;iy++){

    for(int ix=1;ix<N-1;ix++){

      if(ix==1)
        {
            d[ix]=(p[ix][iy]-p[ix+1][iy])*dy-u[ix][iy+1]*(a_n(ix,iy,1))-u[ix][iy-1]*(a_s(ix,iy,1))-u[ix-1][iy]*a_w(ix,iy,1)+(dx*dy/dt)*u_n_1[ix][iy]+fx_n[ix][iy]*dx*dy;
            l[ix]=0; c[ix]=a_e(ix,iy,1); b[ix]=a(ix,iy,1);
            continue;
        }
        if(ix==N-2)
        {
            d[ix]=(p[ix][iy]-p[ix+1][iy])*dy -u[ix][iy+1]*(a_n(ix,iy,1))-u[ix][iy-1]*(a_s(ix,iy,1))-u[ix+1][iy]*a_e(ix,iy,1)+(dx*dy/dt)*u_n_1[ix][iy]+fx_n[ix][iy]*dx*dy;
            l[ix]=a_w(ix,iy,1); b[ix]=a(ix,iy,1);
            continue;
        }
        l[ix]=a_w(ix,iy,1);
        c[ix]=a_e(ix,iy,1);
        b[ix]=a(ix,iy,1);
        d[ix]=(p[ix][iy]-p[ix+1][iy])*dy-u[ix][iy+1]*(a_n(ix,iy,1))-u[ix][iy-1]*(a_s(ix,iy,1))+(dx*dy/dt)*u_n_1[ix][iy]+fx_n[ix][iy]*dx*dy;
                  }
    TDMAn();
    for(int ix=1;ix<N-1;ix++)
        u[ix][iy]=x[ix];

}
//X-Sweep Substitution
/*for(int ix=1;ix<N-1;ix++){

    for(int iy=1;iy<M;iy++){
         if(iy==1)
        {
            d[iy]=(p[ix][iy]-p[ix+1][iy])*dy - u[ix+1][iy]*(a_e(ix,iy,1))-u[ix-1][iy]*(a_w(ix,iy,1))-u[ix][iy-1]*a_s(ix,iy,1);
            l[iy]=0; c[iy]=a_n(ix,iy,1); b[iy]=a(ix,iy,1);
            continue;
        }
        if(iy==M-1)
        {
            d[iy]=(p[ix][iy]-p[ix+1][iy])*dy- u[ix+1][iy]*(a_e(ix,iy,1))-u[ix-1][iy]*(a_w(ix,iy,1))-u[ix][iy+1]*a_n(ix,iy,1);
            l[iy]=a_s(ix,iy,1); b[iy]=a(ix,iy,1);c[iy] = 0;
            continue;
        }
        l[iy]=a_s(ix,iy,1);
        //cout<<l[iy]<<endl;
        c[iy]=a_n(ix,iy,1);
        b[iy]=a(ix,iy,1);
        d[iy]=(p[ix][iy]-p[ix+1][iy])*dy- u[ix+1][iy]*(a_e(ix,iy,1))-u[ix-1][iy]*(a_w(ix,iy,1));
                  }
    TDMAm();
    for(int iy=1;iy<=M-1;iy++)
        {u[ix][iy]=xm[iy];
        //cout<<u[ix][iy]<<endl;
        }
}

*/
}

void gs_p(void){
//Y-Sweep Substitution
for(int iy=1;iy<M;iy++){

   for(int ix=1;ix<N;ix++){
     if(ix==1)
        {
            d[ix]=v[ix][iy-1]*dx-v[ix][iy]*dx+u[ix-1][iy]*dy-u[ix][iy]*dy -p_error[ix][iy-1]*p_s(ix,iy) - p_error[ix][iy+1]*p_n(ix,iy)-p_error[ix-1][iy]*p_w(ix,iy);
            l[ix]=0; c[ix]=p_e(ix,iy); b[ix]=-1.0*(p_n(ix,iy)+p_w(ix,iy)+p_e(ix,iy)+p_s(ix,iy));
            continue;
        }
        if(ix==N-1)
        {
            d[ix]=v[ix][iy-1]*dx-v[ix][iy]*dx+u[ix-1][iy]*dy-u[ix][iy]*dy-p_error[ix][iy-1]*p_s(ix,iy) - p_error[ix][iy+1]*p_n(ix,iy)-p_error[ix+1][iy]*p_e(ix,iy);
            l[ix]=p_w(ix,iy); b[ix]=-1.0*(p_n(ix,iy)+p_w(ix,iy)+p_e(ix,iy)+p_s(ix,iy));
            continue;
        }
        l[ix]=p_w(ix,iy);
        c[ix]=p_e(ix,iy);
        b[ix]=-1.0*(p_n(ix,iy)+p_w(ix,iy)+p_e(ix,iy)+p_s(ix,iy));
        d[ix]=v[ix][iy-1]*dx-v[ix][iy]*dx+u[ix-1][iy]*dy-u[ix][iy]*dy-p_error[ix][iy-1]*p_s(ix,iy) - p_error[ix][iy+1]*p_n(ix,iy);
        }
    TDMApm();
    for(int ix=1;ix<N;ix++)
        p_error[ix][iy]=x[ix];
}
//X-Sweep Substitution
/*
for(int ix=1;ix<N;ix++){

    for(int iy=1;iy<M;iy++){
         if(iy==1)
        {
            d[iy]=v[ix][iy-1]*dx-v[ix][iy]*dx+u[ix-1][iy]*dy-u[ix][iy]*dy -p_error[ix-1][iy]*p_w(ix,iy) - p_error[ix+1][iy]*p_e(ix,iy)-p_error[ix][iy-1]*p_s(ix,iy);
            l[iy]=0; c[iy]=p_n(ix,iy); b[iy]=-1.0*(p_n(ix,iy)+p_w(ix,iy)+p_e(ix,iy)+p_s(ix,iy));
            continue;
        }
        if(iy==M-1)
        {
            d[iy]=v[ix][iy-1]*dx-v[ix][iy]*dx+u[ix-1][iy]*dy-u[ix][iy]*dy -p_error[ix-1][iy]*p_w(ix,iy) - p_error[ix+1][iy]*p_e(ix,iy)-p_error[ix][iy+1]*p_n(ix,iy);
            l[iy]=p_s(ix,iy); b[iy]=-1.0*(p_n(ix,iy)+p_w(ix,iy)+p_e(ix,iy)+p_s(ix,iy));
            continue;
        }
        l[iy]=p_s(ix,iy);
        c[iy]=p_n(ix,iy);
        b[iy]=-1.0*(p_n(ix,iy)+p_w(ix,iy)+p_e(ix,iy)+p_s(ix,iy));
        d[iy]=v[ix][iy-1]*dx-v[ix][iy]*dx+u[ix-1][iy]*dy-u[ix][iy]*dy -p_error[ix-1][iy]*p_w(ix,iy) - p_error[ix+1][iy]*p_e(ix,iy);
                  }
    TDMAm();
    for(int iy=1;iy<M;iy++)
        p_error[ix][iy]=x[iy];
}*/
}


void gs_v(void){

for(int ix=1;ix<N;ix++){

    for(int iy=1;iy<M-1;iy++){

        v[ix][iy]= (dx*1.0*((p[ix][iy]-p[ix][iy+1]))+fy_n[ix][iy]*dx*dy-1.0*(v[ix+1][iy]*(a_e(ix,iy,2)) + v[ix][iy+1]*(a_n(ix,iy,2)) + v[ix-1][iy]*(a_w(ix,iy,2)) + v[ix][iy-1]*(a_s(ix,iy,2)))+v_n_1[ix][iy]*dx*dy/dt)/(a(ix,iy,2));

    }
}

}
//Y-Sweep Elimination
void TDMAn(void){
 int i;
 double t;
    for (i=2;i<=N-2;i++)
    {
        l[i]=l[i]/(b[i-1]);
        b[i]=b[i]-l[i]*c[i-1];
        d[i]=d[i]-l[i]*d[i-1];
    }
    int n=N-2;
    x[n]=d[n]/b[n];
    for (i=n-1;i>=1;i--)
    {
        x[i]=(d[i]-c[i]*x[i+1])/b[i];
        t=x[i];
    }
}
//X-Sweep Elimination
void TDMAm(void){
 int i;
    for (i=2;i<=M-1;i++)
    {
        l[i]=l[i]/(b[i-1]);
        b[i]=b[i]-l[i]*c[i-1];
        d[i]=d[i]-l[i]*d[i-1];
    }
    int n=M-1;
    xm[n]=d[n]/b[n];
    for (i=n-1;i>=1;i--)
    {
        xm[i]=(d[i]-c[i]*xm[i+1])/b[i];
    }
}

void TDMApm(void){
 int i;
    for (i=2;i<=N-1;i++)
    {
        l[i]=l[i]/(b[i-1]);
        b[i]=b[i]-l[i]*c[i-1];
        d[i]=d[i]-l[i]*d[i-1];
    }
    int n=N-1;
    x[n]=d[n]/b[n];
    for (i=n-1;i>=1;i--)
    {
        x[i]=(d[i]-c[i]*x[i+1])/b[i];
    }
}

double UpdateMembranePosition(Lagrange L, int k, int rq)
{
	double sum = 0;
    float factor;
    double xn = 0;
    //ofstream fout398;
      //  fout398.open("Delta.dat");

    if(k == 1)
    {
        for(int i = 1; i < N; i++)
        {
            for(int j = 1; j < M; j++)
            {
                if(rq == 1)
                    sum = sum + u[i][j]*deltah(i*dx - L.xLag, j*dy + 0.5*dy - L.yLag, 1)*dx*dy*0.5*dt;
                if(rq == 2)
                    sum = sum + u[i][j]*deltah(i*dx - L.xLag, j*dy + 0.5*dy - L.yLag, 1)*dx*dy*dt;

            }
        }

        xn = L.xLag + sum;
    }
    double sum1 = 0;
    if(k == 2)
    {
        for(int i = 1; i < N; i++)
        {
            for(int j = 1; j < M; j++)
            {
                if(rq == 1)
                    {
                        sum = sum + v[i][j]*deltah(i*dx + 0.5*dx - L.xLag, j*dy - L.yLag, 1)*dx*dy*0.5*dt;
                        //sum1 = sum1 + deltah(i*dx + 0.5*dx - L.xLag, j*dy - L.yLag, 1);

                    }

                if(rq == 2)
                    sum = sum + v[i][j]*deltah(i*dx + 0.5*dx - L.xLag, j*dy - L.yLag, 1)*dx*dy*dt;
                    //fout398<<v[i][j]<<" "<<deltah(i*dx + 0.5*dx - L.xLag, j*dy - L.yLag, 1)<<" "<<v[i][j]*deltah(i*dx + 0.5*dx - L.xLag, j*dy - L.yLag, 1)*dx*dy*dt<<endl;

            }
        }
       // cout<<"sum1 = "<<sum1*dx*dy<<endl;

        xn = L.yLag + sum;
       // cout<<"xn  =  "<<sum<<endl;
        }

    return  xn;

}


double ForceCalculate(Lagrange L, Lagrange Ln, Lagrange Lp, Lagrange** A, int k)
{
	double Fs = 0, Fb;
    double sum = 0;
    double K;
    double value1, value2, Dx1, Dx2, Dy1, Dy2, ds1, ds2, ds, T1, T2, taux1, taux2, tauy1, tauy2, Fsx, Fsy, F;
	ds1 = sqrt(0.25*(L.xLag - Ln.xLag)*(L.xLag - Ln.xLag) + 0.25*(L.yLag - Ln.yLag)*(L.yLag - Ln.yLag));
	ds2 = sqrt(0.25*(L.xLag - Lp.xLag)*(L.xLag - Lp.xLag) + 0.25*(L.yLag - Lp.yLag)*(L.yLag - Lp.yLag));

    if(L.positionx <= 15 || L.positionx >= NLagMem - 15)
        K = 50000;
    else
        K = Ks;


    /*Dx1 = (Ln.xLag - L.xLag )/dsMem;
    Dx2 = (L.xLag - Lp.xLag )/dsMem;
    Dy1 = (Ln.yLag - L.yLag )/dsMem;
    Dy2 = (L.yLag - Lp.yLag )/dsMem;

        //cout<<value1<<"  "<<L.d_s(Ln)<<endl;
    value1 = sqrt(Dx1*Dx1+Dy1*Dy1);
    value2 = sqrt(Dx2*Dx2+Dy2*Dy2);
    taux1 = Dx1/value1;
    taux2 = Dx2/value1;
    tauy1 = Dy1/value2;
    tauy2 = Dy2/value2;
    T1 = Ks * (abs(value1) - 1);
    T2 = Ks * (abs(value2) - 1);
    Fsx = (T1*taux1-T2*taux2)/dsMem;
    Fsy = (T1*tauy1-T2*tauy2)/dsMem;*/

    int i = L.positionx;
    int j = L.positiony;

    Fsx = K*(L.Zx - L.xLag);
   	Fsy = K*(L.Zy - L.yLag);

    if(k == 1)
    {
           if(i == 0)

              sum = -Kb*(A[L.positiony][1].xLag+A[L.positiony][3].xLag-2*A[L.positiony][2].xLag - (A[L.positiony][0].xLag+A[L.positiony][2].xLag-2*A[L.positiony][1].xLag));
            else if(i == 1)
                sum = -Kb*(A[L.positiony][1].xLag+A[L.positiony][3].xLag-2*A[L.positiony][2].xLag - 2*(A[L.positiony][0].xLag+A[L.positiony][2].xLag-2*A[L.positiony][1].xLag));
            else if(i == NLagMem-1)
                sum = -Kb*(A[L.positiony][i-1].xLag+A[L.positiony][i-3].xLag-2*A[L.positiony][i-2].xLag - (A[L.positiony][i].xLag+A[L.positiony][i-2].xLag-2*A[L.positiony][i-1].xLag));
            else if(i == NLagMem-2)
                sum = -Kb*(A[L.positiony][i].xLag+A[L.positiony][i-2].xLag-2*A[L.positiony][i-1].xLag - 2*(A[L.positiony][i+1].xLag+A[L.positiony][i-1].xLag-2*A[L.positiony][i].xLag));
            else
                sum = -Kb*(A[L.positiony][i+2].xLag - 4*A[L.positiony][i+1].xLag + 6*A[L.positiony][i].xLag - 4*A[L.positiony][i-1].xLag + A[L.positiony][i-2].xLag);
        Fb = sum/(dsMem*dsMem*dsMem*dsMem);
        F = -Fb + Fsx;
        //cout<<"Fbx = "<<Fb<<" "<<"Fsx = "<<Fsx<<endl;
    }

    if(k == 2)
    {
            if(i == 0)
                sum =-Kb*(A[L.positiony][1].yLag+A[L.positiony][3].yLag-2*A[L.positiony][2].yLag - (A[L.positiony][0].yLag+A[L.positiony][2].yLag-2*A[L.positiony][1].yLag));
            else if(i == 1)
                sum = -Kb*(A[L.positiony][1].yLag+A[L.positiony][3].yLag-2*A[L.positiony][2].yLag - 2*(A[L.positiony][0].yLag+A[L.positiony][2].yLag-2*A[L.positiony][1].yLag));
            else if(i == NLagMem-1)
                sum = -Kb*(A[L.positiony][i-1].yLag+A[L.positiony][i-3].yLag-2*A[L.positiony][i-2].yLag - (A[L.positiony][i].yLag+A[L.positiony][i-2].yLag-2*A[L.positiony][i-1].yLag));
            else if(i == NLagMem-2)
                sum = -Kb*(A[L.positiony][i].yLag+A[L.positiony][i-2].yLag-2*A[L.positiony][i-1].yLag - 2*(A[L.positiony][i+1].yLag+A[L.positiony][i-1].yLag-2*A[L.positiony][i].yLag));
            else
             sum = -Kb*(A[L.positiony][i+2].yLag - 4*A[L.positiony][i+1].yLag + 6*A[L.positiony][i].yLag - 4*A[L.positiony][i-1].yLag + A[L.positiony][i-2].yLag);
            Fb = sum/(dsMem*dsMem*dsMem*dsMem);
            F = -Fb + Fsy;
       // cout<<"Fby = "<<Fb<<" "<<"Fsy = "<<Fsy<<endl;
    }
    return F;
}
double EulerianForce(Lagrange**A, int p, int q, int k)
{
	double sum = 0;
	double f;

	if(k == 1)
	{
        for(int i = 0; i < 2; i++)
			for(int j = 1; j < NLagMem - 1; j++)
			{{sum = sum +  A[i][j].Fx_Lag*deltah(p*dx  -  A[i][j].xLag, q*dy + 0.5*dy -  A[i][j].yLag, 1)*dsMem;}}
        	f = sum;
	}

	if(k == 2)
	{
        for(int i = 0; i < 2; i++)
			for(int j = 1; j < NLagMem - 1; j++)
			{{sum = sum +  A[i][j].Fy_Lag*deltah(p*dx + 0.5*dx -  A[i][j].xLag, q*dy  -  A[i][j].yLag, 1)*dsMem;}}

		f = sum;
	}

	return f;
}

double LagrangianVelocity(Lagrange L, int k)
{
	double sum = 0;
	double xn;
	for(int i = 1; i < N - 1; i++)
	{
		for(int j = 1; j < M - 1; j++)
		{
			if(k == 1)
				sum = sum + u[i][j]*deltah((i*dx + 0.5*dx - L.xLag), (j*dy + 0.5*dy - L.yLag), 1);
			if(k == 2)
				sum = sum + v[i][j]*deltah((i*dx + 0.5*dx - L.xLag), (j*dy + 0.5*dy - L.yLag), 1);
		}
	}
	xn = L.yLag + (dx*dy*dt*0.5)*sum;
	return xn;
}
void SIMPLE(int mn, double t)
{
    int b;
    b = length3*(M-1) + 1;
    double f[N+1][M+1];
    double R[b], I[b];
    int r = 0;
    int im = 0;
    ifstream fin1, fin2;
    ofstream fout1;

    fin1.open("Real Part.txt");
    fin2.open("Imaginary Part.txt");

    while(!fin1.eof())
    {
        fin1>>R[r];
        r++;
    }

    while(!fin2.eof())
    {
        fin2>>I[im];
        im++;
    }

     ofstream fout15;
            fout15.open("tolerance.dat");
    if (mn == 1)
        dt = dt/2;

    if (mn == 2)
        dt = dt;
    for(int it=NT;it>0;it--)
    		{

        		cout<<"Iteration number: "<<NT- it<<endl;
               // r = 0;
               //im = 0;
               // fout1.open("Vel Profile.dat");
          		//Boundary Conditions-Velocity
               /* for (int i=(M-1)*(length2 - length3)*0.5;i < (M-1)*(length2 + length3)*0.5;i++)
                {

                    u[0][i] = -(R[r]*sin(2*3.142857*t) + I[im]*cos(3.142857*t));
                    r++;
                    im++;
                    fout1<<u[0][i]<<endl;
                }*/

    			for (int i=1;i<M;i++)
    			{
        			u[0][i] = 0;//Left_Velocity Inlet
        			u[N-1][i]= 0;
        			//u[N-1][i]=u[N-2][i];//Right
        		}
        		/*for (int i=Jin;i<=Jout;i++)
    			{
        			//u[0][i] = 0;//Left_Velocity Inlet
        			u[Iout][i]= u[Iout-1][i];
        			//u[N-1][i]=u[N-2][i];//Right
        		}*/
        		for (int i=1;i<N;i++)
        		{
        			v[i][0]=0.0;//Bottom
        			v[i][M-1]=0.0;//Top
        		}
    			for (int i=0;i<=M-1;i++)
    			{
        			v[0][i]=-v[1][i];//Left
        			v[N][i]=-v[N-1][i];//Right
        		}
    			for (int i=0;i<=N-1;i++)
    			{
        			u[i][0]=-u[i][1];//Bottom
        			u[i][M]=-u[i][M-1];//Top
        		}

        		//Boundary Conditions-Pressure
         		for(int ix=0;ix<=M;ix++)
         		{
            		p[N][ix]=p[N-1][ix];//Right
        		}
        		for(int ix=0;ix<=N;ix++)
        		{
            		p[ix][M]=p[ix][M-1];//Top
        		}
         		for(int ix=0;ix<=M;ix++)
         		{
            		p[0][ix]=p[1][ix];//Left
        		}
        		for(int ix=0;ix<=N;ix++)
        		{
            		p[ix][0]=p[ix][1];//Bottom
        		}
				//p[0][0]=0.0;

				//copying values of U-velocity of previous iteration
        		for(int ix=0;ix<N;ix++)
        		{
            		for(int iy=0;iy<=M;iy++)
            		{
                		un[ix][iy]=u[ix][iy];
            		}
       			}
				//copying values of V-velocity of previous iteration

       			for(int ix=0;ix<=N;ix++)
       			{
        			for(int iy=0;iy<M;iy++)
        			{
            			vn[ix][iy]=v[ix][iy];
        			}
       			}


        		// solving the system of linear equations for u and v momentum equations and pressure correction equation
        		gs_u();
        		gs_v();

        		for(int ix=0;ix<N;ix++)
        		{
            		for(int iy=0;iy<M;iy++)
            		{
                		p_error[ix][iy]=0.0;
            		}
       			}

       			//Solving Pressure Correction Equation
        		gs_p();

          		//updating pressure with under-relaxation factor

            	for(int ix=1;ix<N;ix++)
            	{
                	for(int iy=1;iy<M;iy++)
                	{
                    	p[ix][iy]= p[ix][iy]+ 0.7 * p_error[ix][iy];// here '0.7' is the under relaxation  factor
                   	}
            	}

            	//update velocity values using correction equations

            	//updating u velocity

            	for(int ix=N-2;ix>0;ix--)
            	{
                	for(int iy=M-1;iy>0;iy--)
                	{
                        u[ix][iy]= u[ix][iy]+ (dy*1.0/a(ix,iy,1))*(p_error[ix][iy]-p_error[ix+1][iy]);//updating velocity values without using under relaxation factor
                    }
                }
                //updating v velocity
            	for(int ix=N-1;ix>0;ix--)
            	{
                	for(int iy=M-2;iy>0;iy--)
                	{
                		v[ix][iy]= v[ix][iy]+ (dx*1.0/a(ix,iy,2))*(p_error[ix][iy]-p_error[ix][iy+1]);//updating velocity values without using under relaxation factor
                	}
            	}
            	//Printing values to check if code is working
            	cout<<u[N/2][M/2]<<endl;
            	cout<<u[N][M/2]<<endl;


            	//checking convergence

            	double maximum1=0.0;
            	double maximum2=0.0;
            	double maximum=0.0;

            	for(int ix=0;ix<N-1;ix++)
            	{
                	for(int iy=0;iy<M;iy++)
                	{
    	                if(abs(u[ix][iy]-un[ix][iy])>maximum1)
    	                {
                        	maximum1= abs(u[ix][iy]-un[ix][iy]);//u is velocity in current iteration; un is velocity of previous iteration
                       	}
         	       	}
            	}

            	for(int ix=0;ix<N;ix++)
            	{
                	for(int iy=0;iy<M-1;iy++)
                	{
    		           	if(abs(v[ix][iy]-vn[ix][iy])>maximum1)
    		           	{
                        	maximum1= abs(v[ix][iy]-vn[ix][iy]);
                        }
                	}
            	}

            	if(maximum1<maximum2)
            	{
                	maximum = maximum2;
            	}
            	else if(maximum2<maximum1)
            	{
                	maximum = maximum1;
            	}
            	if(maximum<tolerance)
                break;
    			cout<<"Value of error:"<<maximum<<endl;//Printing value of error after each iteration
    			fout15<<NT-it<<" "<<maximum<<endl;
    		}
    	if (mn == 1)
                dt=2*dt;
}
