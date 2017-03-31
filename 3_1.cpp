//  Graphic;  Vectors; Threads
//  Non-Linear Interconnected Heat exchenger and Mass exchanger Mathematical Model
#include <iostream>
#include <cstring>
#include <math.h>
#include <iomanip>
#include <cstdlib>
#include <GL/glut.h>
#include <vector>
#include <thread>
#include <cstdarg>
#include <stdio.h>
#include <fstream>
#include <time.h>

#define N 100000
#define h 0.5
#define z 5
#define dt 3                    // dt = 3;  //dt = 1; - it's perfect!
#define md 'G'
 
using namespace std;

long double timeBeg = clock();

void TVThread(vector<vector<double> > &TV, double initLayerTV[]);   // Vapor temperature
void TFThread(vector<vector<double> > &TF, double initLayerTF[]);   // Fluid temperature
void CVThread(vector<vector<double> > &CV, double initLayerCV[]);   // Vapor concentration
void CFThread(vector<vector<double> > &CF, double initLayerCF[]);   // Fluid concentration


//void x1Thread(int N,short z, double **TV, double buf1[]);
//void y1Thread(int N,short z, double **y1, double buf2[]);


void output(GLfloat x, GLfloat y, char const *text)
{
    glPushMatrix();
    glTranslatef(x, y, 0);
	GLfloat ficks = 0.1f;
	glScalef(ficks, ficks, 0);
	  
    for( char const *p = text; *p; p++)
    {
        glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, *p);
    }
    glPopMatrix();
}


void displayMe()
{
    // -----Model's heat exchenger parameters------
    double RvT = 0.0191806, RfT = 0.0000777,  a0 = 0.00016966, // a0 = 0.0001152735759,~ ? 0.000207008
    PTV_L = (a0 * 273.15 * dt) / h, PTV_N = 0, PTF = (0.0002291314 * dt) / h,
    initLayerTV[z] = {160, 160, 156, 151.99, 147.99},
    initLayerTF[z] = {120.37, 120.37, 124.38, 128.38, 132.39};

    // -----Model's mass exchenger parameters------
    double RvM = 0.004302, RfM = 0.00001222, E = 0.000000001,
    PCV = (0.07453 * dt) / h, PCF = (0.0002402 * dt) / h,
    initLayerCV[z] = {67.94, 67.94, 69.96, 72.04, 72.04},
    initLayerCF[z] = {6.5, 6.5, 4.613, 2.78, 2.78};


    vector <double> bmp(z, 0);
    vector <vector <double> > TV(N,bmp);
    vector <vector <double> > TF(N,bmp);
    vector <vector <double> > CV(N,bmp);
    vector <vector <double> > CF(N,bmp);

    thread t0(TVThread, ref(TV), initLayerTV);
    thread t1(TFThread, ref(TF), initLayerTF);
    thread t2(CVThread, ref(CV), initLayerCV);
    thread t3(CFThread, ref(CF), initLayerCF);

    t0.join();
    t1.join();
    t2.join();
    t3.join();

/*
    double **TV1 = new double *[N];
    double **y1 = new double *[N];
    double **grphpnts = new double *[N];
	  for(int i = 0; i < N; i++)
    {
        TV[i] = new double [z];
        y1[i] = new double [z];
        grphpnts[i] = new double [z];
    }

    thread t0(x1Thread, N, z, TV, buf1);
    thread t9(y1Thread, N, z, y1, buf2);

    t0.join();
    t9.join();


    for(int i = 0; i < N; i++)
    {
       TV[i][0] = 160;  // 2(1)
       TV[i][1] = 0; TV[i][2] = 0; TV[i][3] = 0;
       TV[i][z-1] = 147.99;
       y1[i][0] = 120.37;  // 2(2)
       y1[i][1] = 0; y1[i][2] = 0; y1[i][3] = 0;
       y1[i][z-1] = 132.39;
    }


    for(int i = 0; i < 5; i++)  // 3(1)
    {
       TV[0][i] = buf1[i];
    }

    for(int i = 0; i < 5; i++) // 3(2)
    {
       y1[0][i] = buf2[i];
    }
*/

    
    // Calculate model
    for(unsigned int i = 1; i < N; i++)   // n: t
    {
       for(short j = 1; j < (z-1); j++)  //  i: z
       {
       		// -----Calculate layer heat exchenger model------
        	PTV_N = (a0*TV[i-1][j+1] * dt) / h;
        	TV[i][j] = -TV[i-1][j] * (PTV_L - 1 - PTV_N + dt*RvT) + (PTV_L * TV[i-1][j-1]) - (PTV_N * TV[i-1][j+1]) + (dt * RvT * TF[i-1][(z-1)-j]);
        	TF[i][j] = -TF[i-1][j] * (PTF - 1 + dt*RfT) + (PTF * TF[i-1][j-1]) + (dt * RfT * TV[i-1][(z-1)-j]);

        	// -----Calculate layer mass exchenger model------
        	CV[i][j] = -CV[i-1][j] * (PTV_L - 1 - PTV_N - dt*RvM) + (PTV_L * CV[i-1][j-1]) -(PTV_N * CV[i-1][j+1]) - (dt * RvM * E * CF[i-1][(z-1)-j]);
			CF[i][j] = -CF[i-1][j] * (PTF - 1 - (dt*RfM*E)) + (PTF * CF[i-1][j-1]) - (dt * RfM * CV[i-1][(z-1)-j]);   	   
       }
    }


	//---------------------------------Switch need function:----------------------------------
	
	GLfloat grphpnts[N][z];

	float cP = 4.358974358, kP = 0.66666, cG = 5.132075471, gS = 0;  // coefficients for scalability
	char const *T0, *T1, *T2;

	if (md == 'P')
	{
		for(unsigned int i = 0; i < N; i++) // 1(0)
		{
		  for(short j = 0; j < z; j++)
		   {
			   grphpnts[i][j] = cP * CV[i][j];
		   }
		   cout << endl;
		}
		gS = 1;
		T0 = "67.94", T1 = "69.96", T2 = "72.04";
	}
	else
	{
	  	if (md == 'G')
	  	{
		    for(unsigned int i = 0; i < N; i++) // 1(0)
		    {
		        for(short j = 0; j < z; j++)
		        {
			         grphpnts[i][j] = cG * CF[i][j];
		        }
		        cout << endl;
		    }
		    gS = 0.1;
		    T0 = "6.5", T1 = "4.613", T2 = "2.78";
    	}
	  	else
	  	{
			cout << "Error input!!!" << endl;
		  	return;
	  	}
	}


	//------------------------------------------------------------------------------------------------
	
	
    glClear(GL_COLOR_BUFFER_BIT);
    glBegin(GL_LINES);

	// Drawing axis coordinates:
	  
    short x0 = 25, xN = 1000, y0 = 5, yN = 680;
	    
    glVertex2f(x0, y0);							// one point of coordinates axis x
    glVertex2f(xN, y0);

	glVertex2f(xN, y0);							// one part point of axis x
    glVertex2f(xN-10, y0+5);
	glVertex2f(xN, y0);							// two part point of axis x
    glVertex2f(xN-10, y0-5);

    glVertex2f(x0, y0);							// one point of coordinates axis y
	glVertex2f(x0, yN);
    
	glVertex2f(x0, yN);							// one part point of axis y
    glVertex2f(x0-4, yN-10);
	glVertex2f(x0, yN);							// two part point of axis y
    glVertex2f(x0+4, yN-10);		


    // Drawing lines on axis y:     ( x = 5.132075471 - coeff scalability;  156 * x = 680 ); 5.132075471*0.66666*grphpnts[0][1]

	glVertex2f( x0+3, kP*grphpnts[0][1] );		// T0
    glVertex2f( x0-3, kP*grphpnts[0][1] );

	glVertex2f( x0+3, grphpnts[0][2]/2 );	    // T1
    glVertex2f( x0-3, grphpnts[0][2]/2 );

	glVertex2f( x0+3, grphpnts[0][3]/3 );		// T2
    glVertex2f( x0-3, grphpnts[0][3]/3 );


	for(short i = 25; i < xN; i += 59)
	{
	    glVertex2f(i, y0+3);					// other lines on axis x
        glVertex2f(i, y0-3);
	}
	
	// Drawing function to lines:
		
    for(short i = 1; i < (z-1); ++i)  // 1(1)
	{
		double xbuf = 25;

		for(unsigned int m = 0; m < N; ++m)
		{
			if (i == 1)										     // for scalability
				glVertex2f(xbuf, grphpnts[m][i]*kP/i);
			else 
				glVertex2f(xbuf, grphpnts[m][i]/i);	

			xbuf += gS;									     // step

			if (i == 1)
				glVertex2f(xbuf, grphpnts[m+1][i]*kP/i);
	        else 
				glVertex2f(xbuf, grphpnts[m+1][i]/i);
		}
    }
		
	
		
	glEnd();

	char yCoord[] = {'C', ',', 'g', '\0'};    // we don't forget symbol for end of string
	char xCoord[] = {'t', '*', '1', '0', '^', '3', ' ', 's', 'e', 'c', '\0'};
	char cCoord[] = {'0', '\0'};
	output(35, 660, yCoord);  // 2(0)
	output(880, 50, xCoord);
	output(10, 5, cCoord);
	
	output(x0+1, (kP*grphpnts[0][1])+12 , T0);
	output(x0+1, (grphpnts[0][2]/2)+12 , T1);
	output(x0+1, (grphpnts[0][3]/3)+12 , T2);
	
	char it[]={"0.625 1.25  1.875 2.5  3.125 3.75 4.375  5  5.625  6.25  6.875 7.5  8.125 8.75 9.375  10"};    // supposably for two functions
	output(x0+(35), y0+10,  it);


	glFlush();


	cout << "Difference time: " << (double)(clock() - timeBeg) / CLOCKS_PER_SEC << endl;

	//----------NNMMTE--------------
	ofstream foutV("Tv.txt"); // 4(1)
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < z; j++)
    	{
      		foutV << TV[i][j] << " | ";
    	}
      	foutV <<"  : " << i*dt << endl;
  	}
  	foutV.close();
	
  	cout << endl << endl;

  	ofstream foutF("Tf.txt");  // 4(2)
  	for(int i = 0; i < N; i++)
  	{
    	for(int j = 0; j < z; j++)
    	{
      		foutF << TF[i][j] << " | ";
    	}
    	foutF <<"  : " << i*dt << endl;
  	}
  	foutF.close();
	
  	cout << endl << endl;

	//----------NNMMME--------------

	ofstream foutCV("Cv.txt"); // 4(1)
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < z; j++)
    	{
      		foutCV << CV[i][j] << " | ";
    	}
      	foutCV <<"  : " << i*dt << endl;
  	}
  	foutCV.close();
	
  	cout << endl << endl;

  	ofstream foutCF("Cf.txt");  // 4(2)
  	for(int i = 0; i < N; i++)
  	{
    	for(int j = 0; j < z; j++)
    	{
      		foutCF << CF[i][j] << " | ";
    	}
    	foutCF <<"  : " << i*dt << endl;
  	}
  	foutCF.close();
	
  	cout << endl << endl;




}



int main(int argc, char** argv)
{
	
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE);
    glutInitWindowSize(800, 400);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Simulation");

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, 1000, 0, 680);

    glutDisplayFunc(displayMe);
    glutMainLoop();

		
    return 0;
}











void TVThread(vector <vector <double> > &TV, double initLayerTV[])
{
    for(int i = 0; i < N; i++)
    {
       TV[i][0] = 160;  // 2(1)
       TV[i][1] = 0; TV[i][2] = 0; TV[i][3] = 0;
       TV[i][z-1] = 147.99;
    }

    for(short i = 0; i < z; i++)  // 3(1)
    {
       TV[0][i] = initLayerTV[i];
    }
    //cout << "1" << endl;
}

void TFThread(vector <vector <double> > &TF, double initLayerTF[])
{
    for(int i = 0; i < N; i++)
    {
       TF[i][0] = 120.37;  // 2(2)
       TF[i][1] = 0; TF[i][2] = 0; TF[i][3] = 0;
       TF[i][z-1] = 132.39;
    }

    for(short i = 0; i < z; i++) // 3(2)
    {
       TF[0][i] = initLayerTF[i];
    }
    //cout << "2" << endl;
}

void CVThread(vector <vector <double> > &CV, double initLayerCV[])
{
    for(int i = 0; i < N; i++)
    {
       CV[i][0] = 67.94;  // 2(1)
       CV[i][1] = 0; CV[i][2] = 0; CV[i][3] = 0;
       CV[i][z-1] = 72.04;
    }

    for(short i = 0; i < z; i++)  // 3(1)
    {
       CV[0][i] = initLayerCV[i];
    }
    //cout << "1" << endl;
}

void CFThread(vector <vector <double> > &CF, double initLayerCF[])
{
    for(int i = 0; i < N; i++)
    {
       CF[i][0] = 6.5;  // 2(2)
       CF[i][1] = 0; CF[i][2] = 0; CF[i][3] = 0;
       CF[i][z-1] = 2.78;
    }

    for(short i = 0; i < z; i++) // 3(2)
    {
       CF[0][i] = initLayerCF[i];
    }
    //cout << "2" << endl;
}



/*

void x1Thread(int N, short z, double **TV, double buf1[])
{
    for(int i = 0; i < N; i++)
    {
       TV[i][0] = 160;  // 2(1)
       TV[i][1] = 0; TV[i][2] = 0; TV[i][3] = 0;
       TV[i][z-1] = 147.99;
    }

    for(short i = 0; i < z; i++)  // 3(1)
    {
       TV[0][i] = buf1[i];
    }
    cout << "1" << endl;
}

void y1Thread(int N, short z, double **y1, double buf2[])
{
    for(int i = 0; i < N; i++)
    {
       y1[i][0] = 120.37;  // 2(2)
       y1[i][1] = 0; y1[i][2] = 0; y1[i][3] = 0;
       y1[i][z-1] = 132.39;
    }

    for(short i = 0; i < z; i++) // 3(2)
    {
       y1[0][i] = buf2[i];
    }
    cout << "2" << endl;
}
*/
