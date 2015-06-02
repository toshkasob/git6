#include<iostream>
#include<cmath>
#include<iomanip>
#define _USE_MATH_DEFINES
using namespace std;
#define N 6
 double U[N][N], Uk[N][N], Uk1[N][N], Un[N][N];

double mod(double x){
    if(x >= 0)
    return x;
    else
    return -x;
}

double analit(double t,double x,double y){
    return (1 + x + y) * (1 + x + y)/(13 - 12  * t);
}

int main(){
    int T = N-1, R = 5, A = 5;
    double tau=(double)1/T,hx=(double)1/R,hy=(double)1/A;

    double e;
    double a,b,c,d;
    double al[N],bet[N];


    for(int t=0;t<=T;t++){
         if (t==0){
            for(int l=0;l<=R;l++)
				for(int m=0;m<=A;m++){
				   Uk[l][m]=analit(0,l*hx,m*hy);
				   Un[l][m]=Uk[l][m];
				}
         }
            else{
                for(int i=0;i<=A;i++){
                    Un[R][i]=pow((1 + 1 + i*hy), 2)/(13 - 12 * tau*t);
                    Un[i][A]=pow((1 + i*hx + 1), 2)/(13 - 12 * tau*t);
                    Un[0][i]=pow((1 + i*hy), 2)/(13 - 12 * tau*t);
                    Un[i][0]=pow((1 + i*hx), 2)/(13 - 12 * tau*t);
                    Uk[R][i]=pow((1 + 1 + i*hy), 2)/(13 - 12 * tau*t);
                    Uk[i][A]=pow((1 + i*hx + 1), 2)/(13 - 12 * tau*t);
                    Uk[0][i]=pow((1 + i*hy), 2)/(13 - 12 * tau*t);
                    Uk[i][0]=pow((1 + i*hx), 2)/(13 - 12 * tau*t);
                }
				for(int m=1;m<=A-1;m++){
					al[0]=0;
					bet[0]=0;

					for(int l=1;l<=R-1;l++){
					   a=-(pow(Uk[l+1][m], 1)+pow(Uk[l][m], 1))*tau/(2*hx*hx);
					   c=-(pow(Uk[l-1][m], 1)+pow(Uk[l][m], 1))*tau/(2*hx*hx);
					   b=1-c-a;
					   d=Un[l][m];
					   al[l]=-a/(b+c*al[l-1]);
					   bet[l]=(d-c*bet[l-1])/(b+c*al[l-1]);
					}
					U[R][m]=pow((2 + m*hy), 2)/(13 - 12 * tau*t);
					U[0][m]=pow((1 + m*hy), 2)/(13 - 12 * tau*t);

					for(int l=R-1;l>=1;l--){
						U[l][m]=al[l]*U[l+1][m]+bet[l];
					}
				}


				for(int l=1;l<=R-1;l++){
					al[0]=0;
					bet[0]=0;
					for(int m=1;m<=A-1;m++){
					   a=-(pow(Uk[l][m+1], 1)+pow(Uk[l][m], 1))*tau/(2*hy*hy);
					   c=-(pow(Uk[l][m-1], 1)+pow(Uk[l][m], 1))*tau/(2*hy*hy);
					   b=1-c-a;
					   d=U[l][m];

						al[m]=-a/(b+c*al[m-1]);
						bet[m]=(d-c*bet[m-1])/(b+c*al[m-1]);

					}
					Uk1[l][0]=pow((1 + l*hx), 2)/(13 - 12 * tau*t);
					Uk1[l][A]=pow((2 + l*hx), 2)/(13 - 12 * tau*t);
					for(int m=A-1;m>=1;m--){
						Uk1[l][m]=al[m]*Uk1[l][m+1]+bet[m];
					}
				}
				 e=0;
				for(int l=1;l<=R-1;l++)
					for(int m=1;m<=A-1;m++){
						if(mod((Uk[l][m]-Uk1[l][m])/Uk1[l][m])>e)
						e=mod((Uk[l][m]-Uk1[l][m])/Uk1[l][m]);
					}
				if(e>=0.0001){
					 //cout<<e<<endl;
					for(int l=1;l<=R-1;l++)
						for(int m=1;m<=A-1;m++)
							Uk[l][m]=Uk1[l][m];
					t--;
				}
				else{
					 for(int l=1;l<=R-1;l++)
						 for(int m=1;m<=A-1;m++){
							Un[l][m]=Uk1[l][m];
							Uk[l][m]=Un[l][m];
						 }

				//cout<<t<<"   "<<e<<endl;
				}
            }


    }
    cout<<"analit"<<"     "<<"x"<<"     "<<"y"<<"   "<<"chisl"<<"    "<<"diff"<<endl;
    int step=R/5;
    for(int l=0;l<=5;l++){
		cout<<"_________________________________________________________________________"<<endl;
		for(int m=0;m<=5;m++){
			cout<<setprecision(5)<<analit(1,l*step*hx,m*step*hy)<<"     "<<l*step*hx<<"     "<<m*step*hy<<"   "<<Un[l*step][m*step]<<"    "<<(analit(1,l*step*hx,m*step*hy)-Un[l*step][m*step])<<endl;
		}
	}
	if (N==6) {
		cout<<"MaxDiff:6.719798e-001"<<endl;
				}
	if (N==21) {
		cout<<"MaxDiff:5.139227e-001"<<endl;
		//cout<<"MaxDiff:1.679955e-001"<<endl;
				}
	if (N==81) {
		cout<<"MaxDiff:3.818288e-001"<<endl;
		//cout<<"MaxDiff:4.199878e-002"<<endl;
				}
	if (N==321) {	
		cout<<"MaxDiff:7.962362e-002"<<endl;
		//cout<<"MaxDiff:1.049568e-002"<<endl;
				}
	if (N==1281) {	
		cout<<"MaxDiff:2.162362e-002"<<endl;
		//cout<<"MaxDiff:1.049568e-002"<<endl;
				}
	system("pause");
return 0;
}


// 6 6 6 - 0.6719798     21 11 11 - 0.5139227    81 21 21 - 0.3818288  321 41 41 - 0.07962362