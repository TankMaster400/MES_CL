#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>      
using namespace std;

double ff2(double x, double y)
{
    return 5 * pow(x, 2) * pow(y, 2) + 3 * x * y + 6;
}
double ff1(double x)
{
    return 5 * pow(x, 2) + 3 * x + 6;
}

double gauss(int num, int dim)
{
    double* X = new double[num];
    double* W = new double[num];

    switch (num)
    {
    case 2:
        X[0] = -sqrt(3.0) / 3.0;
        X[1] = sqrt(3.0) / 3.0;
        W[0] = 1.0;
        W[1] = 1.0;
        break;
    case 3:
        X[0] = -sqrt(15.0) / 5.0;
        X[1] = 0.0;
        X[2] = sqrt(15.0) / 5.0;
        W[0] = 5.0 / 9.0;
        W[1] = 8.0 / 9.0;
        W[2] = 5.0 / 9.0;
        break;
    case 4:
        X[0] = -sqrt(525 + 70 * sqrt(30)) / 35.0;
        X[1] = -sqrt(525 - 70 * sqrt(30)) / 35.0;
        X[2] = sqrt(525 - 70 * sqrt(30)) / 35.0;
        X[3] = sqrt(525 + 70 * sqrt(30)) / 35.0;
        W[0] = (18.0 - sqrt(30.0)) / 36.0;
        W[1] = (18.0 + sqrt(30.0)) / 36.0;
        W[2] = (18.0 + sqrt(30.0)) / 36.0;
        W[3] = (18.0 - sqrt(30.0)) / 36.0;
        break;
    }


    double wart = 0;

    switch (dim)
    {
    case 1:

        for (int i = 0; i < num; i++)
        {
            wart += W[i] * ff1(X[i]);
        }
        return wart;
    case 2:
        for (int i = 0; i < num; i++)
        {
            for (int j = 0; j < num; j++)
            {
                wart += W[i]* W[j] * ff2(X[i], X[j]);
            }
        }
        return wart;
    }

}

double dN1e(double n)
{
    return -0.25 * (1 - n);
}

double dN2e(double n)
{
    return 0.25 * (1 - n);
}

double dN3e(double n)
{
    return 0.25 * (1 + n);
}

double dN4e(double n)
{
    return -0.25 * (1 + n);
}

double dN1n(double e)
{
    return -0.25 * (1 - e);
}

double dN2n(double e)
{
    return -0.25 * (1 + e);
}

double dN3n(double e)
{
    return 0.25 * (1 + e);
}

double dN4n(double e)
{
    return 0.25 * (1 - e);
}


struct Global_data
{
    int SimulationTime;
    int SimulationStepTime;
    int Conductivity;
    int Alfa;
    int Tot;
    int InitialTemp;
    int Density;
    int SpecificHeat;
    Global_data() {};
    Global_data(int ST, int SST, int C, int A, int T, int IT, int D, int SH): SimulationTime(ST),SimulationStepTime(SST),Conductivity(C), Alfa(A),Tot(T),InitialTemp(IT),Density(D),SpecificHeat(SH){}
};
struct node
{
    double x,y;
    int t ;
    node();
};
struct element
{
    int ID[1][4];
    element();
};

struct J
{
    double** Tab;
    J(int  n)
    {
        Tab = new double* [n];
   

        for (int i = 0; i < n; i++)
        {
            Tab[i] = new double[2];
        }
        int i = sqrt(n);
     
        switch (i)
        {
        case 2:
            Tab[0][0] = -sqrt(3.0) / 3.0;
            Tab[0][1] = -sqrt(3.0) / 3.0;
            Tab[1][0] = -sqrt(3.0) / 3.0;
            Tab[1][1] = sqrt(3.0) / 3.0;
            Tab[2][0] = sqrt(3.0) / 3.0;
            Tab[2][1] = -sqrt(3.0) / 3.0;
            Tab[3][0] = sqrt(3.0) / 3.0;
            Tab[3][1] = sqrt(3.0) / 3.0;
          
            break;
        case 3:
            
            Tab[0][0] = -sqrt(15.0) / 5.0;
            Tab[0][1] = -sqrt(15.0) / 5.0;
            Tab[1][0] = 0;
            Tab[1][1] = -sqrt(15.0) / 5.0;
            Tab[2][0] = sqrt(15.0) / 5.0;
            Tab[2][1] = -sqrt(15.0) / 5.0;
            Tab[3][0] = -sqrt(15.0) / 5.0;
            Tab[3][1] = 0;
            Tab[4][0] = 0;
            Tab[4][1] = 0;
            Tab[5][0] = sqrt(15.0) / 5.0;
            Tab[5][1] = 0;
            Tab[6][0] = -sqrt(15.0) / 5.0;
            Tab[6][1] = sqrt(15.0) / 5.0;
            Tab[7][0] = 0;
            Tab[7][1] = sqrt(15.0) / 5.0;
            Tab[8][0] = sqrt(15.0) / 5.0;
            Tab[8][1] = sqrt(15.0) / 5.0;
     
            break;
        }

    }
};

struct element_uni
{
    double ** Tab;
    double ** Tab2;
    int N;

    element_uni(int n, struct J wart_p)
    {
        Tab = new double*[n];
        Tab2 = new double*[n];
        N = n;

        for (int i = 0; i < n; i++)
        {
            Tab[i] = new double[4];
        }

        for (int i = 0; i < n; i++)
        {
                    
                    Tab[i][0] = dN1e(wart_p.Tab[i][0]);

                    Tab[i][1] = dN2e(wart_p.Tab[i][0]);

                    Tab[i][2] = dN3e(wart_p.Tab[i][0]);

                    Tab[i][3] = dN4e(wart_p.Tab[i][0]);
        }

        for (int i = 0; i < n; i++)
        {
            Tab2[i] = new double[4];
        }

        for (int i = 0; i < n; i++)
        {
            Tab2[i][0] = dN1n(wart_p.Tab[i][1]);

            Tab2[i][1] = dN2n(wart_p.Tab[i][1]);

            Tab2[i][2] = dN3n(wart_p.Tab[i][1]);

            Tab2[i][3] = dN4n(wart_p.Tab[i][1]);
        }
  }
    
};

struct matrixH
{
    double** dNdx;
    double** dNdy;
    double detJ [4];
    double tk[4][2][2];
    double H[4][4][4];
    matrixH(double ID_xy[2][4], struct element_uni el)
    {
       for(int i = 0; i < 4 ; i++)
       { 
        tk[i][0][0] = el.Tab[i][0] * ID_xy[0][0] + el.Tab[i][1] * ID_xy[0][1] + el.Tab[i][2] * ID_xy[0][2] + el.Tab[i][3] * ID_xy[0][3];
        tk[i][0][1] = el.Tab[i][0] * ID_xy[1][0] + el.Tab[i][1] * ID_xy[1][1] + el.Tab[i][2] * ID_xy[1][2] + el.Tab[i][3] * ID_xy[1][3];
        tk[i][1][0] = el.Tab2[i][0] * ID_xy[0][0] + el.Tab2[i][1] * ID_xy[0][1] + el.Tab2[i][2] * ID_xy[0][2] + el.Tab2[i][3] * ID_xy[0][3];
        tk[i][1][1] = el.Tab2[i][0] * ID_xy[1][0] + el.Tab2[i][1] * ID_xy[1][1] + el.Tab2[i][2] * ID_xy[1][2] + el.Tab2[i][3] * ID_xy[1][3];

        detJ[i] = 1 / (tk[i][1][1] * tk[i][0][0] - tk[i][0][1] * tk[i][1][1]);
       }
        dNdx = new double* [el.N];
        dNdy = new double* [el.N];
 

        for (int i = 0; i < el.N; i++)
        {
            dNdx[i] = new double[4];
        }

        for (int i = 0; i < el.N; i++)
        {
            dNdy[i] = new double[4];
        }

        for (int i = 0; i < el.N; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                dNdx[i][j] = tk[i][0][0] * detJ[i] * el.Tab[i][j] + tk[i][1][0] * detJ[i] * el.Tab2[i][j];
            }
        }

        for (int i = 0; i < el.N; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                dNdy[i][j] = tk[i][0][1] * detJ[i] * el.Tab[i][j] + tk[i][1][1] * detJ[i] * el.Tab2[i][j];
            }
        }
       
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                for (int g = 0; g < 4; g++)
                {
                    H[i][j][g] = 30 * (dNdx[i][j] * dNdx[i][g] + dNdy[i][j] * dNdy[i][g]) * 1 / detJ[i];
                }
            
            }
        }
    }
};

struct grid
{
    int Nn;
    int En;
    
     node * Tnode ;
     element * Tele ;
     grid(int N, int E) : Nn(N), En(E) {  };
};



int main()
{
    ifstream readfile("Test1_4_4.txt");
    double valgd[8];

    for (int i = 0; i < 8; i++)
    {
        readfile.ignore(std::numeric_limits<std::streamsize>::max(), ' ') && readfile >> valgd[i];
    }
    
    Global_data GD(valgd[0], valgd[1], valgd[2], valgd[3], valgd[4], valgd[5], valgd[6],  valgd[7]);
    
    cout << "Global data" << endl;
    cout << GD.SimulationTime << endl;
    cout << GD.SimulationStepTime << endl;
    cout << GD.Conductivity << endl;
    cout << GD.Alfa<< endl;
    cout << GD.Tot << endl;
    cout << GD.InitialTemp << endl;
    cout << GD.Density << endl;
    cout << GD.SpecificHeat << endl;

    int n, e;

    readfile.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
    readfile.ignore(std::numeric_limits<std::streamsize>::max(), ' ') && readfile >> n;
    readfile.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
    readfile.ignore(std::numeric_limits<std::streamsize>::max(), ' ') && readfile >> e;

    grid grid1(n, e);
    grid1.Tnode = (node*)malloc(sizeof(node) * grid1.Nn);
    grid1.Tele = (element*) malloc(sizeof(element) * grid1.En);
    
     double test, test2;
    for (int i = 0; i < grid1.Nn; i++)
    {
        readfile.ignore(std::numeric_limits<std::streamsize>::max(), ',' ) && readfile >> test;
        readfile.ignore(std::numeric_limits<std::streamsize>::max(), ',') && readfile >> test2;

        grid1.Tnode[i].x = test;
        grid1.Tnode[i].y = test2;
        grid1.Tnode[i].t = 1;
    }

    readfile.ignore(std::numeric_limits<std::streamsize>::max(), ',');

    int e1;

    for (int i = 0; i < grid1.En; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            readfile.ignore(std::numeric_limits<std::streamsize>::max(), ',') && readfile >> e1;
            grid1.Tele[i].ID[0][j] = e1;
        }

    }

    cout << "GRID" << endl;
    cout << grid1.En << endl;
    cout << grid1.Nn << endl;
    for (int i = 0; i < grid1.Nn; i++)
    {
        
        cout << i +1 <<"   " << std::setprecision(9) << grid1.Tnode[i].x << "                     " << grid1.Tnode[i].y << endl;
       
    }
    for (int i = 0; i < grid1.En; i++)
    {
        cout << i + 1;
        for (int j = 0; j < 4 ;j++)
        {
          cout<< "   " << grid1.Tele[i].ID[0][j] ;
        }
        cout << endl;
    }

    cout << gauss(2, 1) << endl;
    cout << gauss(2, 2) << endl;
    cout << gauss(3, 1) << endl;
    cout << gauss(3, 2) << endl;

    int N= 4;

    J Jac(N);

    element_uni el(N, Jac);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            cout << "   " << el.Tab[i][j];
        }
        cout << endl;

    }
    cout << endl;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            cout << "   " << el.Tab2[i][j];
        }
        cout << endl;

    }
    cout << endl;
    double ID_xy[2][4] = { {0,0.025,0.025, 0},{0,0,0.025,0.025} };
    matrixH  H(ID_xy, el);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            cout << "   " << H.dNdx[i][j];
        }
        cout << endl;

    }
    cout << endl;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            cout << "   " << H.dNdy[i][j];
        }
        cout << endl;

    }

    cout << endl;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            for (int g = 0; g < 4; g++)
            {
                cout << "   " << H.H[i][j][g];
            }
            cout << endl;
        }
        cout << endl;
    }
    free(grid1.Tnode);
    free(grid1.Tele);
    return 0;
}
