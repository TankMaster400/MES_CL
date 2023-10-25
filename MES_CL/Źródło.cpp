#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>      
using namespace std;
//test afafda
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

double dN1(double n)
{
    return -0.25 * (1 - n);
}

double dN2(double n)
{
    return 0.25 * (1 - n);
}

double dN3(double n)
{
    return 0.25 * (1 + n);
}

double dN4(double n)
{
    return -0.25 * (1 + n);
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
struct element_uni
{
    double * Tab[4];

    element_uni(int n)
    {
        for (int i = 0; i < 4; i++)
        {
            Tab[i] = new double[n];
        }

        for (int i = 0; i < 4; i++)
        {
            switch (i)
            {
            case 0:

                for (int j = 0; j < n; j++)
                {

                    Tab[i][j] = dN1(i);
                }
                break;
            case 1:

                for (int j = 0; j < n; j++)
                {

                    Tab[i][j] = dN2(i);
                }
                break;
            case 2:

                for (int j = 0; j < n; j++)
                {

                    Tab[i][j] = dN3(i);
                }
                break;
            case 3:

                for (int j = 0; j < n; j++)
                {

                    Tab[i][j] = dN4(i);
                }
                break;
        
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

    element_uni el(4);

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            cout << "   " << el.Tab[i][j];
        }
        cout << endl;
    }


    free(grid1.Tnode);
    free(grid1.Tele);
    return 0;
}
