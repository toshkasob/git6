//sorry, it's not working because :ya zaputalsya: :D

// _______________________
// what need to do?
// try change. what is done?
// it's work. changing is highlighted by green in commits
// _______________________
// translit=vozmozhnie oshibki ili ih ispravlenie
// clear English = real comment to understand how program work
// _______________________
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
// #include <ifstream>
// #include <ofstream>
#define PI 3.14159265
#define EPS 0.0001
using namespace std;

int H, T, F;
float h, t, f;					// step by r, t, fi
float mu, Ct, delta, tT, fF;			//tT=t*T:0<t<tT		//fF=f*F:0<f<fF		//delta is used to standart tT=1
enum WherePoint
{ inside, outside };
float analit(float x, float f, float y)
{
	WherePoint flag;
	flag = (x < 0 || x > 1 || f < 0 || f > PI / 1.9 || y < 0 || y > 1) ? outside : inside;
	switch (flag)
	{
	case inside:
		// cout << "belowLC" << endl;
		return (x * cos(f) / sqrt(5 - 4 * y));
		// solution 
		break;
	default:
		cout << "INCORRECT! point x=" << x << "  t=" << y << "fi=" << f << 
			"  is not in [0,1]×[0,1]" << endl;
		throw 0;
		return -10;				// all solution >=
		break;
	}
}





int main(){
	cout.setf(ios::showpos);
	char restart = 'n';
	mu=2; Ct=5;
	H = 5;
	T = 5;
	F = 5;
	delta=1./4.;
	tT=Ct*mu/(2*mu+2*2)-delta;
	fF=PI/2;
	float maxEps = 0;
  start:
	t = tT / T;
	h = 1. / H;
	f = fF / F;
	float px, py;
	/* { cout << "enter the point x="; cin >> px; cout << " y="; cin >> py;
	   cout << endl; float mindSolution = analit(px, py); // sled if
	   (mindSolution == -10) { cout << "error. аналитическое
	   решение не может быть найдено... конец
	   программы" << endl; return 1; } else { cout << "mindSolution=" 
	   << mindSolution << endl; }

	   } */

/*	vector < float > mindSolutionAtT(H + 1);		//budet vichislyat'sya pered vivodom
	for (int i = 0; i <= H; i++)
	{
		try
		{
			mindSolutionAtT[i] = analit(i * h, 1);
		}
		catch(int f)
		{
			cout <<
				"error. your mindsolution cannot calculated. you need change X in your programm";
			return 0;
		}
	}
	cout << endl;
*/
	vector <float> alfa(H+1);
	vector <float> betta(H+1);
//	vector <float> cl(H+1);
//vector <float> dl(H+1);
	
	vector < vector  < float > > un(H + 1);
	vector < vector  < float > > uk(H + 1);
	vector < vector  < float > > U(H + 1);
	vector < vector  < float > > U1(H + 1);
// NB: i*h=xL, q*t=tN, j*f=fi
	for (int i = 0; i <= H; i++){
		un[i].resize(F + 1);
		uk[i].resize(F + 1);
		U[i].resize(F + 1);
		U1[i].resize(F + 1);
	}
	

//_filling of the original data
	for (int q=0; q<T; q++){

		if (!q){
			for (int i = 0; i <= H; i++) {
				for (int j = 0; j <= F; j++){
					un[i][j]=analit(i*h,j*f,0);
					uk[i][j]=un[i][j];
				}
			}
		} else {
			for (int m=0;m<=F && m<=H;m++){
				un[0][m]=0;
				un[m][0]=analit(1,m*f,q*t);
				un[H][m]=analit(m*h,0,q*t);
				un[m][F]=0;
				uk[0][m]=0;
				uk[m][0]=analit(1,m*f,q*t);
				uk[H][m]=analit(m*h,0,q*t);
				uk[m][F]=0;
			}
//__the first half-step	
			for (int j=1;j<=F-1;j++){
				alfa[0]=0;
				betta[0]=0;
				for (int i=1;i<=H-1;i++){
					float a=-h*(i+0.5)*(pow(uk[i+1][j],mu)+pow(uk[i][j],mu))*t/(2*(i*h)*h*h);
					float c=-h*(i-0.5)*(pow(uk[i][j],mu)+pow(uk[i-1][j],mu))*t/(2*(i*h)*h*h);
					float b=1-a-c;
					float d=un[i][j];
					alfa[i]=-a/(b+c*alfa[i-1]);
					betta[i]=(d-c*betta[i-1])/(b+c*alfa[i-1]);
				}
				U[0][j]=0;
				U[H][j]=analit(1,j*f,(q+1)*t);
				for(int i=H-1;i>=1;i--){
					U[i][j]=alfa[i]*U[i+1][j]+betta[i];
				}
			}
//___the second half-step	
			for (int i=1;i<=H-1;i++){
				alfa[0]=0;
				betta[0]=0;
				for (int j=1;j<=F-1;j++){
					float a=-(pow(uk[i][j+1],mu)+pow(uk[i][j],mu))*t/(2*(i*h)*(i*h)*f*f);
					float c=-(pow(uk[i][j],mu)+pow(uk[i][j-1],mu))*t/(2*(i*h)*(i*h)*f*f);
					float b=1-a-c;
					float d=U[i][j];
					alfa[j]=-a/(b+c*alfa[j-1]);
					betta[j]=(d-c*betta[j-1])/(b+c*alfa[j-1]);
				}
				U1[i][F]=0;
				U1[i][0]=analit(i*h,0,(q+1)*t);
				for(int j=F-1;j>=1;j--){
					U1[i][j]=alfa[i]*U[i+1][j]+betta[i];
				}
			}
//____check(29)
			maxEps=0;
			for(int i=1;i<=H-1;i++){
				for(int j=1;j<=F-1;j++){
					if(fabs( (U1[i][j]-uk[i][j]) / U1[i][j])>maxEps)
						maxEps=fabs( (U1[i][j]-uk[i][j]) / U1[i][j]);
				}
			}
			if (maxEps>=EPS){
				for(int i=1;i<=H-1;i++){
					for(int j=1;j<=F-1;j++){
						uk[i][j]=U1[i][j];
					}
				}
				q--;
			} else {
				for(int i=1;i<=H-1;i++){
					for(int j=1;j<=F-1;j++){
						un[i][j]=U1[i][j];
						uk[i][j]=un[i][j];
					}
				}
			}
//_____

				
		}
		

	}
//______ cout final table:
	int columnWidth = 10;
	float maxDiff = 0.;
	cout << endl;
	cout.precision(4);
	cout << "_____" << setw(columnWidth*2) << "sled. " << "_____"<< endl ;
	for(int i=-1;i<=H;i++){
		for (int j=-1;j<=F;j++){
			if (i==-1){
				if (j==-1)
					cout << setw(columnWidth )<< "r\fi ";
				else 
					cout << setw(columnWidth/2) << "fi=" << j*f;
			} else {
				if (j==-1)
					cout << setw(columnWidth/2) << "r=" << i*h;
				else
					cout << setw(columnWidth) << analit(i*h,j*f,1);
			}
		}
		cout<<endl;
	}


	cout << "_____" << setw(columnWidth*2) << "chislennoe reshenie." << "_____" << endl;
	for(int i=-1;i<=H;i++){
		for (int j=-1;j<=F;j++){
			if (i==-1){
				if (j==-1)
					cout << setw(columnWidth) << "r\fi ";
				else 
					cout << setw(columnWidth/2) << "fi=" << j*f;
			} else {
				if (j==-1)
					cout << setw(columnWidth/2) << "r=" << i*h;
				else
					cout << setw(columnWidth) << un[i][j];
			}
		}
		cout<<endl;
	}


	cout << "_____" << setw(columnWidth*2) << "difference module. " << "_____" << endl;
	for(int i=-1;i<=H;i++){
		for (int j=-1;j<=F;j++){
			if (i==-1){
				if (j==-1)
					cout << setw(columnWidth) << "r\fi ";
				else 
					cout << setw(columnWidth/2) << "fi=" << j*f;
			} else {
				if (j==-1)
					cout << setw(columnWidth/2) << "r=" << i*h;
				else {
					float diff=fabs(un[i][j]-analit(i*h,j*f,1));
					if (maxDiff < diff){
						maxDiff = diff;
					}
					cout << setw(columnWidth) << diff;
				}
			}
		}
		cout<<endl;
	}


	/*cout << setw(columnWidth / 2) << "x. ";
	cout << setw(columnWidth) << "sled. " ;//- analiticheskoe reshenie"// 
	cout << setw(columnWidth) << "chislennoe reshenie. ";
	cout << setw(columnWidth) << "difference module. " << endl;
	for (int i = 0; i <= H; i++)
	{
		if (i % (int)pow(2, flagDouble) == 0)
		{
			// cout.width(columnWidth);
			cout.precision(4);
			cout << setw(columnWidth / 2) << i * h;
			cout << setw(columnWidth) << mindSolutionAtT[i];
			cout << setw(columnWidth) << u[i][T];
			float diff = abs(u[i][T] - mindSolutionAtT[i]);
			cout << setw(columnWidth) << diff;
			if (maxDiff < diff)
			{
				maxDiff = diff;
			}
			cout << endl;
		}
	}*/
	cout << "maxDiff=" << maxDiff << endl;
	system("pause");
/*	// end:
	// cout << "do you want recalc in another point? ";
	cout << "do you want double step?   ";
	cin >> restart;
	if (restart == 'y')
	{
		T *= 2;
		H *= 2;
		flagDouble++;
		goto start;
	}
*/
	return 0;
}