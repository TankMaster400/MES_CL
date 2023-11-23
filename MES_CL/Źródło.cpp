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
    Global_data(int ST, int SST, int C, int A, int T, int IT, int D, int SH) : SimulationTime(ST), SimulationStepTime(SST), Conductivity(C), Alfa(A), Tot(T), InitialTemp(IT), Density(D), SpecificHeat(SH) {}
};

struct node
{
    double x,y;
    double BC;
    node();
};
struct element
{
    int ID[1][4];
    double H[4][4];
    double Hbc[4][4];

 
    element();
   
};

struct GaussIntegration
{
    double** Tab_pc;
    double* Tab_w;

    GaussIntegration(int  n)
    {
        Tab_pc = new double* [n];
        Tab_w = new double [n];

        for (int i = 0; i < n; i++)
        {
            Tab_pc[i] = new double[2];
        }

            Tab_w = new double[n];
        

        int i = sqrt(n);
     
        switch (i)
        {
        case 2:
            Tab_pc[0][0] = -sqrt(3.0) / 3.0;
            Tab_pc[0][1] = -sqrt(3.0) / 3.0;
            Tab_pc[1][0] = -sqrt(3.0) / 3.0;
            Tab_pc[1][1] = sqrt(3.0) / 3.0;
            Tab_pc[2][0] = sqrt(3.0) / 3.0;
            Tab_pc[2][1] = -sqrt(3.0) / 3.0;
            Tab_pc[3][0] = sqrt(3.0) / 3.0;
            Tab_pc[3][1] = sqrt(3.0) / 3.0;
            Tab_w[0] = 1;
            Tab_w[1] = 1;
            Tab_w[2] = 1;
            Tab_w[3] = 1;
          
            break;
        case 3:
            Tab_pc[0][0] = -sqrt(15.0) / 5.0;
            Tab_pc[0][1] = -sqrt(15.0) / 5.0;
            Tab_pc[1][0] = 0;
            Tab_pc[1][1] = -sqrt(15.0) / 5.0;
            Tab_pc[2][0] = sqrt(15.0) / 5.0;
            Tab_pc[2][1] = -sqrt(15.0) / 5.0;
            Tab_pc[3][0] = -sqrt(15.0) / 5.0;
            Tab_pc[3][1] = 0;
            Tab_pc[4][0] = 0;
            Tab_pc[4][1] = 0;
            Tab_pc[5][0] = sqrt(15.0) / 5.0;
            Tab_pc[5][1] = 0;
            Tab_pc[6][0] = -sqrt(15.0) / 5.0;
            Tab_pc[6][1] = sqrt(15.0) / 5.0;
            Tab_pc[7][0] = 0;
            Tab_pc[7][1] = sqrt(15.0) / 5.0;
            Tab_pc[8][0] = sqrt(15.0) / 5.0;
            Tab_pc[8][1] = sqrt(15.0) / 5.0;
            break;
        case 4:
     
            break;
        }

    }
};

struct element_uni
{
    double ** dNdKsi;
    double ** dNdEta;
    int N;
    double S[4];
    double Sur[2][4];

    element_uni(int n, struct GaussIntegration wart_p)
    {
        dNdKsi = new double*[n];
        dNdEta = new double*[n];
        N = n;

        for (int i = 0; i < n; i++)
        {
            dNdKsi[i] = new double[4];
        }

        for (int i = 0; i < n; i++)
        {
                    
            dNdKsi[i][0] = dN1e(wart_p.Tab_pc[i][0]);

            dNdKsi[i][1] = dN2e(wart_p.Tab_pc[i][0]);

            dNdKsi[i][2] = dN3e(wart_p.Tab_pc[i][0]);

            dNdKsi[i][3] = dN4e(wart_p.Tab_pc[i][0]);
        }

        for (int i = 0; i < n; i++)
        {
            dNdEta[i] = new double[4];
        }

        for (int i = 0; i < n; i++)
        {
            dNdEta[i][0] = dN1n(wart_p.Tab_pc[i][1]);

            dNdEta[i][1] = dN2n(wart_p.Tab_pc[i][1]);

            dNdEta[i][2] = dN3n(wart_p.Tab_pc[i][1]);

            dNdEta[i][3] = dN4n(wart_p.Tab_pc[i][1]);
        }
  }
    // tu  cos?

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
    //ifstream readfile("Test1_4_4.txt");
    ifstream readfile("Test2_4_4_MixGrid.txt");
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
        grid1.Tnode[i].BC= 0;
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

    int bc;
    readfile.ignore(std::numeric_limits<std::streamsize>::max(), 'C'); 
    readfile >> bc;
    grid1.Tnode[bc].BC = 1;
    while (readfile.ignore(std::numeric_limits<std::streamsize>::max(), ',') && readfile >> bc)
    {
        grid1.Tnode[bc].BC = 1;
    }
            
            

    cout << "GRID" << endl;
    cout << grid1.En << endl;
    cout << grid1.Nn << endl;
    for (int i = 0; i < grid1.Nn; i++)
    {
        
        cout << i +1 <<"   " << std::setprecision(9) << grid1.Tnode[i].x << "                     " << grid1.Tnode[i].y << "                     " << grid1.Tnode[i].BC<<  endl;
       
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

    int N = 4; //liczba punktów calkowania
    GaussIntegration Jac(N);
    element_uni el(N, Jac);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            cout << "   " << el.dNdKsi[i][j];
        }
        cout << endl;

    }
    cout << endl;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            cout << "   " << el.dNdEta[i][j];
        }
        cout << endl;

    }

    cout << endl;
    double ID_xy[2][4]; // = { {0,0.025,0.025, 0},{0,0,0.025,0.025} };
    double tk[4][2][2];
    double detJ[4];
    double** dNdx;
    double** dNdy;

    dNdx = new double* [N];
    dNdy = new double* [N];

    for (int i = 0; i < N; i++)
    {
        dNdx[i] = new double[4];
    }

    for (int i = 0; i < N; i++)
    {
        dNdy[i] = new double[4];
    }

    for (int i = 0; i < grid1.En; i++)
    {
        cout << i << endl;
        for (int j = 0; j < 4; j++)
        {
            ID_xy[0][j] = grid1.Tnode[grid1.Tele[i].ID[0][j] - 1].x;
            ID_xy[1][j] = grid1.Tnode[grid1.Tele[i].ID[0][j] - 1].y;
        }
        cout << endl;
        for (int j = 0; j < 4; j++)
        {
            tk[j][0][0] = el.dNdKsi[j][0] * ID_xy[0][0] + el.dNdKsi[j][1] * ID_xy[0][1] + el.dNdKsi[j][2] * ID_xy[0][2] + el.dNdKsi[j][3] * ID_xy[0][3];
            tk[j][0][1] = el.dNdKsi[j][0] * ID_xy[1][0] + el.dNdKsi[j][1] * ID_xy[1][1] + el.dNdKsi[j][2] * ID_xy[1][2] + el.dNdKsi[j][3] * ID_xy[1][3];
            tk[j][1][0] = el.dNdEta[j][0] * ID_xy[0][0] + el.dNdEta[j][1] * ID_xy[0][1] + el.dNdEta[j][2] * ID_xy[0][2] + el.dNdEta[j][3] * ID_xy[0][3];
            tk[j][1][1] = el.dNdEta[j][0] * ID_xy[1][0] + el.dNdEta[j][1] * ID_xy[1][1] + el.dNdEta[j][2] * ID_xy[1][2] + el.dNdEta[j][3] * ID_xy[1][3];
            cout << tk[j][0][0] << "   " << tk[j][0][1] << endl;
            cout << tk[j][1][0] << "   " << tk[j][1][1] << endl;
            detJ[j] = (tk[j][1][1] * tk[j][0][0] - tk[j][0][1] * tk[j][1][0]);
            cout << detJ[j] << endl;
            cout << endl;

            for (int g = 0; g < 4; g++)
            {
                dNdx[j][g] = tk[j][0][0] * (1 / detJ[j]) * el.dNdKsi[j][g] + tk[j][1][0] * (1 / detJ[j]) * el.dNdEta[j][g];
            }
            for (int g = 0; g < 4; g++)
            {
                dNdy[j][g] = tk[j][0][1] * (1 / detJ[j]) * el.dNdKsi[j][g] + tk[j][1][1] * (1 / detJ[j]) * el.dNdEta[j][g];
            }
        }
     
        for (int j = 0; j < 4; j++)
        {
            cout << " dNdx " << j << endl;
            for (int g = 0; g < 4; g++)
            {
                cout << dNdx[j][g] <<"  ";
            }
      
            cout << endl;
            cout << " dNdy " << j << endl;
            for (int g = 0; g < 4; g++)
            {
                cout << dNdy[j][g] << "  ";
            }
            cout << endl;
          }
        cout << endl;
        for (int j = 0; j < 4; j++)
        {
            for (int g = 0; g < 4; g++)
            {

                grid1.Tele[i].H[j][g] = 0;
            }
          
        }
        for (int c = 0; c < 4; c++)
        {
            for (int j = 0; j < 4; j++)
            {
                for (int g = 0; g < 4; g++)
                {

                    grid1.Tele[i].H[j][g] += GD.Conductivity * (dNdx[c][j] * dNdx[c][g] + dNdy[c][j] * dNdy[c][g]) * detJ[j];
                   
                }
              
            }
        }
            cout << endl;
            //Dziwne wyniki
            for (int j = 0; j < 4; j++)
            {
                for (int g = 0; g < 4; g++)
                {

                    cout << grid1.Tele[i].H[j][g] << "   ";

                }
                cout << endl;
            }
            

    }
    

    cout << endl;
    return 0;
}
