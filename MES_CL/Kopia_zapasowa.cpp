//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <string>
//#include <iomanip>      
//using namespace std;
//
//double ff2(double x, double y)
//{
//    return 5 * pow(x, 2) * pow(y, 2) + 3 * x * y + 6;
//}
//double ff1(double x)
//{
//    return 5 * pow(x, 2) + 3 * x + 6;
//}
//
//double N1(double ksi, double eta)
//{
//    return 0.25 * (1 - ksi) * (1 - eta);
//}
//double N2(double ksi, double eta)
//{
//    return 0.25 * (1 + ksi) * (1 - eta);
//}
//double N3(double ksi, double eta)
//{
//    return 0.25 * (1 + ksi) * (1 + eta);
//}
//double N4(double ksi, double eta)
//{
//    return 0.25 * (1 - ksi) * (1 + eta);
//}
//
//double dN1e(double n)
//{
//    return -0.25 * (1 - n);
//}
//double dN2e(double n)
//{
//    return 0.25 * (1 - n);
//}
//double dN3e(double n)
//{
//    return 0.25 * (1 + n);
//}
//double dN4e(double n)
//{
//    return -0.25 * (1 + n);
//}
//
//double dN1n(double ksi)
//{
//    return -0.25 * (1 - ksi);
//}
//double dN2n(double ksi)
//{
//    return -0.25 * (1 + ksi);
//}
//double dN3n(double ksi)
//{
//    return 0.25 * (1 + ksi);
//}
//double dN4n(double ksi)
//{
//    return 0.25 * (1 - ksi);
//}
//
//struct Global_data
//{
//    int SimulationTime;
//    int SimulationStepTime;
//    int Conductivity;
//    int Alfa;
//    int Tot;
//    int InitialTemp;
//    int Density;
//    int SpecificHeat;
//    Global_data() {};
//    Global_data(int ST, int SST, int C, int A, int T, int IT, int D, int SH) : SimulationTime(ST), SimulationStepTime(SST), Conductivity(C), Alfa(A), Tot(T), InitialTemp(IT), Density(D), SpecificHeat(SH) {}
//};
//
//struct node
//{
//    double x, y;
//    double BC;
//    node();
//};
//
//struct element
//{
//    int ID[4];
//    double H[4][4];
//    double Hbc[4][4][4];
//    double P[4];
//    double C[4][4];
//    element();
//};
//
//struct GaussIntegration
//{
//    double** Tab_pc;
//    double* Tab_w;
//
//    GaussIntegration(int  n)
//    {
//        Tab_pc = new double* [n];
//        Tab_w = new double[n];
//
//        for (int i = 0; i < n; i++)
//        {
//            Tab_pc[i] = new double[2];
//        }
//
//        int i = sqrt(n);
//        Tab_w = new double[i];
//        switch (i) // 1 - eta || 0 - ksi (ksi, eta) hehe
//        {
//        case 2:
//            Tab_pc[0][0] = -sqrt(3.0) / 3.0;
//            Tab_pc[0][1] = -sqrt(3.0) / 3.0;
//            Tab_pc[1][0] = sqrt(3.0) / 3.0;
//            Tab_pc[1][1] = -sqrt(3.0) / 3.0;
//
//            Tab_pc[2][0] = sqrt(3.0) / 3.0;
//            Tab_pc[2][1] = sqrt(3.0) / 3.0;
//            Tab_pc[3][0] = -sqrt(3.0) / 3.0;
//            Tab_pc[3][1] = sqrt(3.0) / 3.0;
//
//            Tab_w[0] = 1;
//            Tab_w[1] = 1;
//
//
//            break;
//        case 3:
//            Tab_pc[0][0] = -sqrt(15.0) / 5.0;
//            Tab_pc[0][1] = -sqrt(15.0) / 5.0;
//            Tab_pc[1][0] = 0;
//            Tab_pc[1][1] = -sqrt(15.0) / 5.0;
//            Tab_pc[2][0] = sqrt(15.0) / 5.0;
//            Tab_pc[2][1] = -sqrt(15.0) / 5.0;
//
//            Tab_pc[3][0] = -sqrt(15.0) / 5.0;
//            Tab_pc[3][1] = 0;
//            Tab_pc[4][0] = 0;
//            Tab_pc[4][1] = 0;
//            Tab_pc[5][0] = sqrt(15.0) / 5.0;
//            Tab_pc[5][1] = 0;
//
//            Tab_pc[6][0] = -sqrt(15.0) / 5.0;
//            Tab_pc[6][1] = sqrt(15.0) / 5.0;
//            Tab_pc[7][0] = 0;
//            Tab_pc[7][1] = sqrt(15.0) / 5.0;
//            Tab_pc[8][0] = sqrt(15.0) / 5.0;
//            Tab_pc[8][1] = sqrt(15.0) / 5.0;
//
//            Tab_w[0] = 5. / 9.;
//            Tab_w[1] = 8. / 9.;
//            Tab_w[2] = 5. / 9.;
//
//            break;
//        case 4:
//
//            break;
//        }
//
//    }
//};
//
//void HB_PC_wart(int n, double*** Surp)
//{
//    switch (n)
//    {
//    case 2:
//        Surp[0][0][0] = -sqrt(3.0) / 3.0;
//        Surp[0][0][1] = -1.;
//        Surp[0][1][0] = sqrt(3.0) / 3.0;
//        Surp[0][1][1] = -1.;
//
//        Surp[1][0][0] = 1;
//        Surp[1][0][1] = -sqrt(3.0) / 3.0;
//        Surp[1][1][0] = 1.;
//        Surp[1][1][1] = sqrt(3.0) / 3.0;
//
//        Surp[2][0][0] = sqrt(3.0) / 3.0;
//        Surp[2][0][1] = 1.;
//        Surp[2][1][0] = -sqrt(3.0) / 3.0;
//        Surp[2][1][1] = 1.;
//
//        Surp[3][0][0] = -1.;
//        Surp[3][0][1] = sqrt(3.0) / 3.0;
//        Surp[3][1][0] = -1.;
//        Surp[3][1][1] = -sqrt(3.0) / 3.0;
//
//        break;
//    case 3:
//        Surp[0][0][0] = -sqrt(15.0) / 5.0;
//        Surp[0][0][1] = -1.;
//        Surp[0][1][0] = 0;
//        Surp[0][1][1] = -1.;
//        Surp[0][1][0] = sqrt(15.0) / 5.0;
//        Surp[0][1][1] = -1.;
//
//        Surp[1][0][0] = 1;
//        Surp[1][0][1] = -sqrt(15.0) / 5.0;
//        Surp[1][1][0] = 1.;
//        Surp[1][1][1] = 0;
//        Surp[1][1][0] = 1.;
//        Surp[1][1][1] = sqrt(15.0) / 5.0;
//
//        Surp[2][0][0] = -sqrt(15.0) / 5.0;
//        Surp[2][0][1] = 1.;
//        Surp[2][1][0] = 0;
//        Surp[2][1][1] = 1.;
//        Surp[2][0][0] = sqrt(15.0) / 5.0;
//        Surp[2][0][1] = 1.;
//        Surp[2][1][0] = -sqrt(15.0) / 5.0;
//        Surp[2][1][1] = 1.;
//
//        Surp[3][0][0] = -1.;
//        Surp[3][0][1] = sqrt(15.0) / 5.0;
//        Surp[3][1][0] = -1.;
//        Surp[3][1][1] = 0;
//        Surp[3][1][0] = -1.;
//        Surp[3][1][1] = -sqrt(15.0) / 5.0;
//        break;
//    case 4:
//
//        break;
//    }
//};
//
//struct element_uni //tutaj
//{
//    double** dNdKsi;
//    double** dNdEta;
//    double** Tab_N;
//    int N;
//    double S[4];
//    double*** Surp;
//
//    element_uni(int n, struct GaussIntegration wart_p)
//    {
//        dNdKsi = new double* [n];
//        dNdEta = new double* [n];
//        Tab_N = new double* [n];
//        int npc = sqrt(n);
//
//        Surp = new double** [4];
//
//        for (int i = 0; i < 4; i++)
//        {
//            Surp[i] = new double* [npc];
//        }
//        for (int i = 0; i < 4; i++)
//        {
//            for (int j = 0; j < npc; j++)
//            {
//                Surp[i][j] = new double[2];
//            }
//        }
//        for (int i = 0; i < n; i++)
//        {
//            Tab_N[i] = new double[4];
//        }
//
//        N = n;
//
//        for (int i = 0; i < n; i++)
//        {
//            dNdKsi[i] = new double[4];
//        }
//
//        for (int i = 0; i < n; i++)
//        {
//
//            dNdKsi[i][0] = dN1e(wart_p.Tab_pc[i][1]);
//
//            dNdKsi[i][1] = dN2e(wart_p.Tab_pc[i][1]);
//
//            dNdKsi[i][2] = dN3e(wart_p.Tab_pc[i][1]);
//
//            dNdKsi[i][3] = dN4e(wart_p.Tab_pc[i][1]);
//        }
//
//        for (int i = 0; i < n; i++)
//        {
//            dNdEta[i] = new double[4];
//        }
//
//        for (int i = 0; i < n; i++)
//        {
//            dNdEta[i][0] = dN1n(wart_p.Tab_pc[i][0]);
//
//            dNdEta[i][1] = dN2n(wart_p.Tab_pc[i][0]);
//
//            dNdEta[i][2] = dN3n(wart_p.Tab_pc[i][0]);
//
//            dNdEta[i][3] = dN4n(wart_p.Tab_pc[i][0]);
//        }
//
//        for (int i = 0; i < n; i++)
//        {
//            Tab_N[i][0] = N1(wart_p.Tab_pc[i][0], wart_p.Tab_pc[i][1]);
//            Tab_N[i][1] = N2(wart_p.Tab_pc[i][0], wart_p.Tab_pc[i][1]);
//            Tab_N[i][2] = N3(wart_p.Tab_pc[i][0], wart_p.Tab_pc[i][1]);
//            Tab_N[i][3] = N4(wart_p.Tab_pc[i][0], wart_p.Tab_pc[i][1]);
//        }
//
//        HB_PC_wart(npc, Surp);
//    }
//
//};
//
//struct grid
//{
//    int Nn;
//    int En;
//    node* Tnode;
//    element* Tele;
//
//    grid(int N, int E)
//    {
//        Nn = N;
//        En = E;
//        Tnode = (node*)malloc(sizeof(node) * Nn);
//        Tele = (element*)malloc(sizeof(element) * En);
//    }
//};
//
//int main()
//{
//    // ifstream readfile("Test1_4_4.txt");
//    ifstream readfile("Test2_4_4_MixGrid.txt");
//    double valgd[8];
//
//    for (int i = 0; i < 8; i++)
//    {
//        readfile.ignore(std::numeric_limits<std::streamsize>::max(), ' ') && readfile >> valgd[i];
//    }
//
//    Global_data GD(valgd[0], valgd[1], valgd[2], valgd[3], valgd[4], valgd[5], valgd[6], valgd[7]);
//
//    cout << "Global data" << endl;
//    cout << GD.SimulationTime << endl;
//    cout << GD.SimulationStepTime << endl;
//    cout << GD.Conductivity << endl;
//    cout << GD.Alfa << endl;
//    cout << GD.Tot << endl;
//    cout << GD.InitialTemp << endl;
//    cout << GD.Density << endl;
//    cout << GD.SpecificHeat << endl;
//
//    int n, e;
//
//    readfile.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
//    readfile.ignore(std::numeric_limits<std::streamsize>::max(), ' ') && readfile >> n;
//    readfile.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
//    readfile.ignore(std::numeric_limits<std::streamsize>::max(), ' ') && readfile >> e;
//
//    grid grid1(n, e);
//
//    double test, test2;
//    for (int i = 0; i < grid1.Nn; i++)
//    {
//        readfile.ignore(std::numeric_limits<std::streamsize>::max(), ',') && readfile >> test;
//        readfile.ignore(std::numeric_limits<std::streamsize>::max(), ',') && readfile >> test2;
//
//        grid1.Tnode[i].x = test;
//        grid1.Tnode[i].y = test2;
//        grid1.Tnode[i].BC = 0;
//    }
//
//    readfile.ignore(std::numeric_limits<std::streamsize>::max(), ',');
//
//    int e1;
//
//    for (int i = 0; i < grid1.En; i++)
//    {
//        for (int j = 0; j < 4; j++)
//        {
//            readfile.ignore(std::numeric_limits<std::streamsize>::max(), ',') && readfile >> e1;
//            grid1.Tele[i].ID[j] = e1;
//        }
//    }
//
//    int bc;
//    readfile.ignore(std::numeric_limits<std::streamsize>::max(), 'C');
//    readfile >> bc;
//    grid1.Tnode[bc - 1].BC = 1;
//    while (readfile.ignore(std::numeric_limits<std::streamsize>::max(), ',') && readfile >> bc)
//    {
//        grid1.Tnode[bc - 1].BC = 1;
//    }
//
//    cout << "GRID" << endl;
//    cout << grid1.En << endl;
//    cout << grid1.Nn << endl;
//
//    for (int i = 0; i < grid1.Nn; i++)
//    {
//        cout << i + 1 << "   " << std::setprecision(9) << grid1.Tnode[i].x << "                     " << grid1.Tnode[i].y << "                     " << grid1.Tnode[i].BC << endl;
//    }
//    for (int i = 0; i < grid1.En; i++)
//    {
//        cout << i + 1;
//        for (int j = 0; j < 4; j++)
//        {
//            cout << "   " << grid1.Tele[i].ID[j];
//        }
//        cout << endl;
//    }
//    //liczba punktów calkowania
//    int N = 4;
//    //  int N = 9;
//    GaussIntegration GausI(N);
//
//    element_uni el(N, GausI);
//
//    //Wypisanie dNdKsi oraz dNdEta
//
//   /* for (int i = 0; i < N; i++)
//    {
//        for (int j = 0; j < 4; j++)
//        {
//            cout << "   " << el.dNdKsi[i][j];
//        }
//        cout << endl;
//
//    }
//    cout << endl;
//    for (int i = 0; i < N; i++)
//    {
//        for (int j = 0; j < 4; j++)
//        {
//            cout << "   " << el.dNdEta[i][j];
//        }
//        cout << endl;
//
//    }
//    cout << endl;
//    for (int i = 0; i < N; i++)
//    {
//        for (int j = 0; j < 4; j++)
//        {
//            cout << "   " << el.Tab_N[i][j];
//        }
//        cout << endl;
//
//    }*/
//    cout << endl;
//
//    double ID_xy[2][4];  // = { {0,0.025,0.025, 0},{0,0,0.025,0.025} };
//      double*** tk;
//        double* detJ;
//    double** dNdx;
//    double** dNdy;
//    detJ = new double[N];
//    dNdx = new double* [N];
//    dNdy = new double* [N];
//    tk = new double** [N];
//
//    for (int i = 0; i < N; i++)
//    {
//        dNdx[i] = new double[4];
//        dNdy[i] = new double[4];
//        tk[i] = new double* [2];
//    }
//
//    for (int i = 0; i < N; i++)
//    {
//        for (int j = 0; j < 2; j++)
//        {
//            tk[i][j] = new double[2];
//        }
//    }
//
//    double** agr_tab;
//    double** agr_tab_C;
//    double* agr_P;
//    double* agr_C;
//    agr_tab = new double* [grid1.Nn];
//    agr_tab_C = new double* [grid1.Nn];
//    agr_P = new double[grid1.Nn];
//    agr_C = new double[grid1.Nn];
//    for (int j = 0; j < grid1.Nn; j++)
//    {
//        agr_tab[j] = new double[grid1.Nn];
//        agr_tab_C[j] = new double[grid1.Nn];
//    }
//    for (int j = 0; j < grid1.Nn; j++)
//    {
//        for (int g = 0; g < grid1.Nn; g++)
//        {
//            agr_tab[j][g] = 0;
//            agr_tab_C[j][g] = 0;
//        }
//        agr_P[j] = 0;
//        agr_C[j] = 0;
//    }
//
//    //Pêtla elementów
//    for (int i = 0; i < grid1.En; i++)
//    {
//        cout << "Element: " << i << endl;
//        cout << endl;
//
//        //Pobranie punktów dla 
//        for (int j = 0; j < 4; j++)
//        {
//            ID_xy[0][j] = grid1.Tnode[grid1.Tele[i].ID[j] - 1].x;
//            ID_xy[1][j] = grid1.Tnode[grid1.Tele[i].ID[j] - 1].y;
//        }
//
//        //Obliczanie Jakobianu 
//        for (int j = 0; j < N; j++)
//        {
//            tk[j][0][0] = el.dNdKsi[j][0] * ID_xy[0][0] + el.dNdKsi[j][1] * ID_xy[0][1] + el.dNdKsi[j][2] * ID_xy[0][2] + el.dNdKsi[j][3] * ID_xy[0][3];
//            tk[j][0][1] = el.dNdKsi[j][0] * ID_xy[1][0] + el.dNdKsi[j][1] * ID_xy[1][1] + el.dNdKsi[j][2] * ID_xy[1][2] + el.dNdKsi[j][3] * ID_xy[1][3];
//            tk[j][1][0] = el.dNdEta[j][0] * ID_xy[0][0] + el.dNdEta[j][1] * ID_xy[0][1] + el.dNdEta[j][2] * ID_xy[0][2] + el.dNdEta[j][3] * ID_xy[0][3];
//            tk[j][1][1] = el.dNdEta[j][0] * ID_xy[1][0] + el.dNdEta[j][1] * ID_xy[1][1] + el.dNdEta[j][2] * ID_xy[1][2] + el.dNdEta[j][3] * ID_xy[1][3];
//            cout << tk[j][0][0] << "   " << tk[j][0][1] << endl;
//            cout << tk[j][1][0] << "   " << tk[j][1][1] << endl;
//            detJ[j] = (tk[j][1][1] * tk[j][0][0] - tk[j][0][1] * tk[j][1][0]);
//            cout << detJ[j] << endl;
//            cout << endl;
//
//            for (int g = 0; g < N; g++)
//            {
//                dNdx[j][g] = tk[j][1][1] * (1. / detJ[j]) * el.dNdKsi[j][g] + -tk[j][0][1] * (1. / detJ[j]) * el.dNdEta[j][g];
//            }
//            for (int g = 0; g < N; g++)
//            {
//                dNdy[j][g] = -tk[j][1][0] * (1. / detJ[j]) * el.dNdKsi[j][g] + tk[j][0][0] * (1. / detJ[j]) * el.dNdEta[j][g];
//            }
//        }
//        double temp;
//        if (N == 4)
//        {
//            temp = detJ[3];
//            detJ[3] = detJ[2];
//            detJ[2] = temp;
//
//            for (int g = 0; g < 4; g++)
//            {
//                temp = dNdx[3][g];
//                dNdx[3][g] = dNdx[2][g];
//                dNdx[2][g] = temp;
//
//                temp = dNdy[3][g];
//                dNdy[3][g] = dNdy[2][g];
//                dNdy[2][g] = temp;
//            }
//        }
//        //Wypisanie dN/dx oraz dN/dy 
//      /*  for (int j = 0; j < 4; j++)
//        {
//            cout << " dNdx " << j << endl;
//            for (int g = 0; g < 4; g++)
//            {
//                cout << dNdx[j][g] << "  ";
//            }
//
//            cout << endl;
//            cout << " dNdy " << j << endl;
//            for (int g = 0; g < 4; g++)
//            {
//                cout << dNdy[j][g] << "  ";
//            }
//            cout << endl;
//        }
//        cout << endl;*/
//
//        for (int j = 0; j < 4; j++)
//        {
//            for (int g = 0; g < 4; g++)
//            {
//                grid1.Tele[i].H[j][g] = 0;
//            }
//        }
//
//        int w1 = 0, w2 = 0;
//        for (int c = 0; c < N; c++)
//        {
//            if (c == 2) { w1 = 0; }
//            for (int j = 0; j < 4; j++)
//            {
//                for (int g = 0; g < 4; g++)
//                {
//                    grid1.Tele[i].H[j][g] += (GD.Conductivity * (dNdx[c][j] * dNdx[c][g] + dNdy[c][j] * dNdy[c][g]) * detJ[c]) * GausI.Tab_w[w1] * GausI.Tab_w[w2];
//                }
//            }
//            w1++;
//            if (c == 1) { w2++; }
//        }
//
//        cout << "H:" << endl;
//        for (int j = 0; j < 4; j++)
//        {
//            for (int g = 0; g < 4; g++)
//            {
//                cout << grid1.Tele[i].H[j][g] << "   ";
//            }
//            cout << endl;
//        }
//
//        double*** tabN;
//        tabN = new double** [4];
//
//        for (int i = 0; i < 4; i++)
//        {
//            tabN[i] = new double* [sqrt(N)];
//        }
//
//        for (int i = 0; i < 4; i++)
//        {
//            for (int j = 0; j < sqrt(N); j++)
//            {
//                tabN[i][j] = new double[4];
//            }
//        }
//
//        for (int j = 0; j < 4; j++)
//        {
//            for (int g = 0; g < sqrt(N); g++)
//            {
//                tabN[j][g][0] = N1(el.Surp[j][g][0], el.Surp[j][g][1]);
//                tabN[j][g][1] = N2(el.Surp[j][g][0], el.Surp[j][g][1]);
//                tabN[j][g][2] = N3(el.Surp[j][g][0], el.Surp[j][g][1]);
//                tabN[j][g][3] = N4(el.Surp[j][g][0], el.Surp[j][g][1]);
//            }
//            cout << endl;
//        }
//        for (int a = 0; a < 4; a++)
//        {
//            for (int j = 0; j < 4; j++)
//            {
//                for (int g = 0; g < 4; g++)
//                {
//                    grid1.Tele[i].Hbc[a][j][g] = 0;
//                }
//            }
//        }
//        double L;
//
//        for (int a = 0; a < 4; a++)
//        {
//            for (int nn = 0; nn < sqrt(N); nn++)
//            {
//                for (int g = 0; g < 4; g++)
//                {
//                    L = pow((grid1.Tnode[grid1.Tele[i].ID[a % 4] - 1].x - grid1.Tnode[grid1.Tele[i].ID[(a + 1) % 4] - 1].x), 2) + pow((grid1.Tnode[grid1.Tele[i].ID[a % 4] - 1].y - grid1.Tnode[grid1.Tele[i].ID[((a % 4) + 1) % 4] - 1].y), 2);
//                    L = sqrt(L) / 2;
//
//                    for (int c = 0; c < 4; c++)
//                    {
//                        grid1.Tele[i].Hbc[a][g][c] += 300 * GausI.Tab_w[nn] * (tabN[a][nn][g] * tabN[a][nn][c]) * L;
//                    }
//                }
//            }
//        }
//
//        cout << "Hbc:" << endl;
//
//        for (int a = 0; a < 4; a++)
//        {
//            for (int j = 0; j < 4; j++)
//            {
//                for (int g = 0; g < 4; g++)
//                {
//                    cout << grid1.Tele[i].Hbc[a][j][g] << "   ";
//                }
//                cout << endl;
//            }
//            cout << endl;
//        }
//
//        //     cout << endl;
//         //    cout << "P:" << endl;
//        double P[4] = { 0.0,0.0,0.0,0.0 };
//
//        for (int k = 0; k < 4; k++)
//        {
//            for (int c = 0; c < sqrt(N); c++)
//            {
//                for (int g = 0; g < 4; g++)
//                {
//                    L = pow((grid1.Tnode[grid1.Tele[i].ID[k % 4] - 1].x - grid1.Tnode[grid1.Tele[i].ID[((k % 4) + 1) % 4] - 1].x), 2) + pow((grid1.Tnode[grid1.Tele[i].ID[k % 4] - 1].y - grid1.Tnode[grid1.Tele[i].ID[((k % 4) + 1) % 4] - 1].y), 2);
//                    L = sqrt(L) / 2;
//
//                    if (grid1.Tnode[grid1.Tele[i].ID[k % 4] - 1].BC == 1 && grid1.Tnode[grid1.Tele[i].ID[((k % 4) + 1) % 4] - 1].BC == 1)
//                    {
//                        P[g] += 300 * (GausI.Tab_w[c] * tabN[k][c][g] * GD.Tot) * L;
//                    }
//                }
//            }
//        }
//        //cout << " P: " << endl;
//        //for (int a = 0; a < 4; a++)
//        //{
//        //    cout << P[a] << "   ";
//        //}
//        //cout << endl;
//        cout << "H + Hbc" << endl;
//
//        if (grid1.Tnode[grid1.Tele[i].ID[0] - 1].BC == 1 && grid1.Tnode[grid1.Tele[i].ID[1] - 1].BC == 1)
//        {
//            for (int j = 0; j < 4; j++)
//            {
//                for (int g = 0; g < 4; g++)
//                {
//                    grid1.Tele[i].H[j][g] += grid1.Tele[i].Hbc[0][j][g];
//                }
//            }
//        }
//        if (grid1.Tnode[grid1.Tele[i].ID[1] - 1].BC == 1 && grid1.Tnode[grid1.Tele[i].ID[2] - 1].BC == 1)
//        {
//            for (int j = 0; j < 4; j++)
//            {
//                for (int g = 0; g < 4; g++)
//                {
//                    grid1.Tele[i].H[j][g] += grid1.Tele[i].Hbc[1][j][g];
//                }
//            }
//        }
//        if (grid1.Tnode[grid1.Tele[i].ID[2] - 1].BC == 1 && grid1.Tnode[grid1.Tele[i].ID[3] - 1].BC == 1)
//        {
//            for (int j = 0; j < 4; j++)
//            {
//                for (int g = 0; g < 4; g++)
//                {
//                    grid1.Tele[i].H[j][g] += grid1.Tele[i].Hbc[2][j][g];
//                }
//            }
//        }
//        if (grid1.Tnode[grid1.Tele[i].ID[3] - 1].BC == 1 && grid1.Tnode[grid1.Tele[i].ID[0] - 1].BC == 1)
//        {
//            for (int j = 0; j < 4; j++)
//            {
//                for (int g = 0; g < 4; g++)
//                {
//                    grid1.Tele[i].H[j][g] += grid1.Tele[i].Hbc[3][j][g];
//                }
//            }
//        }
//
//        for (int j = 0; j < 4; j++)
//        {
//            for (int g = 0; g < 4; g++)
//            {
//                cout << grid1.Tele[i].H[j][g] << "   ";
//            }
//            cout << endl;
//        }
//        for (int j = 0; j < 4; j++)
//        {
//            for (int g = 0; g < 4; g++)
//            {
//                grid1.Tele[i].C[j][g] = 0;
//            }
//
//        }
//        temp = detJ[2];
//        detJ[2] = detJ[3];
//        detJ[3] = temp;
//        w1 = 0;
//        w2 = 0;
//        cout << "C - " << endl;
//
//        for (int c = 0; c < N; c++)
//        {
//            if (c == 2) { w1 = 0; }
//            for (int j = 0; j < 4; j++)
//            {
//                for (int g = 0; g < 4; g++)
//                {
//                    grid1.Tele[i].C[j][g] += GD.Density * GD.SpecificHeat * (el.Tab_N[c][j] * el.Tab_N[c][g]) * detJ[c] * GausI.Tab_w[w1] * GausI.Tab_w[w2];
//                }
//            }
//            w1++;
//            if (c == 1) { w2++; }
//        }
//
//        for (int j = 0; j < 4; j++)
//        {
//            for (int g = 0; g < 4; g++)
//            {
//
//                cout << grid1.Tele[i].C[j][g] << "  ";
//            }
//            cout << endl;
//        }
//        cout << "C - " << endl;
//
//        for (int j = 0; j < 4; j++)
//        {
//            for (int g = 0; g < 4; g++)
//            {
//                agr_tab_C[grid1.Tele[i].ID[j] - 1][grid1.Tele[i].ID[g] - 1] += grid1.Tele[i].C[j][g];
//            }
//        }
//
//        for (int j = 0; j < 4; j++)
//        {
//            for (int g = 0; g < 4; g++)
//            {
//                agr_tab[grid1.Tele[i].ID[j] - 1][grid1.Tele[i].ID[g] - 1] += grid1.Tele[i].H[j][g];
//            }
//        }
//
//        for (int g = 0; g < 4; g++)
//        {
//            agr_P[grid1.Tele[i].ID[g] - 1] += P[g];
//
//        }
//
//        for (int j = 0; j < 4; j++)
//        {
//            for (int g = 0; g < 4; g++)
//            {
//                agr_C[grid1.Tele[i].ID[g] - 1] += grid1.Tele[i].C[j][g];
//            }
//        }
//        /*  for (int g = 0; g < grid1.Nn; g++)
//          {
//              cout << agr_P[g] << "      ";
//          }*/
//    }
//    cout << endl;
//    cout << std::setprecision(6);
//    for (int j = 0; j < grid1.Nn; j++)
//    {
//        for (int g = 0; g < grid1.Nn; g++)
//        {
//            cout << agr_tab[j][g] + agr_tab_C[j][g] / GD.SimulationStepTime << "      ";
//        }
//        cout << endl;
//    }
//
//    double* agr_PC = new double[grid1.Nn];
//
//    cout << endl;
//    double* tab_wynik = new double[grid1.Nn];
//    double* Y = new double[grid1.Nn];
//    double** L = new double* [grid1.Nn];
//    double** U = new double* [grid1.Nn];
//    for (int i = 0; i < grid1.Nn; i++)
//    {
//        L[i] = new double[grid1.Nn];
//    }
//    for (int i = 0; i < grid1.Nn; i++)
//    {
//        U[i] = new double[grid1.Nn];
//    }
//
//    double Asr = (double)GD.Conductivity / (double)(GD.SpecificHeat * GD.Density);
//
//    double dl = pow((grid1.Tnode[grid1.Tele[2].ID[0] - 1].x - grid1.Tnode[grid1.Tele[2].ID[1] - 1].x), 2) + pow((grid1.Tnode[grid1.Tele[2].ID[0] - 1].y - grid1.Tnode[grid1.Tele[2].ID[1] - 1].y), 2);
//    double dtime = dl / (Asr * 0.5);
//    //dtime = dl / (Asr);
//    //dtime = 67;
//    GD.SimulationStepTime = dtime;
//
//    for (int g = 0; g < grid1.Nn; g++)
//    {
//        agr_PC[g] = agr_P[g] + GD.InitialTemp * agr_C[g] / GD.SimulationStepTime;
//
//    }
//
//    for (int g = 0; g < grid1.Nn; g++)
//    {
//        cout << agr_PC[g] << "      ";
//    }
//
//    for (double p = GD.SimulationStepTime; p < GD.SimulationTime + 1; p += GD.SimulationStepTime)
//    {
//        for (int i = 0; i < grid1.Nn; i++)
//        {
//            for (int j = 0; j < grid1.Nn; j++)
//            {
//                L[i][j] = 0;
//            }
//        }
//        for (int i = 0; i < grid1.Nn; i++)
//        {
//            for (int j = 0; j < grid1.Nn; j++)
//            {
//                U[i][j] = 0;
//            }
//        }
//        for (int i = 0; i < grid1.Nn; i++)
//        {
//            L[i][i] = 1;
//        }
//
//        for (int i = 0; i < grid1.Nn; i++)
//        {
//            for (int j = 0; j < grid1.Nn; j++)
//            {
//                double sum = 0;
//                for (int k = 0; k <= i - 1; k++)
//                {
//                    sum += L[i][k] * U[k][j];
//                }
//
//                if (i <= j)
//                {
//                    U[i][j] = (agr_tab[i][j] + (agr_tab_C[i][j] / GD.SimulationStepTime)) - sum;
//                }
//                else
//                {
//                    L[i][j] = (1. / U[j][j]) * ((agr_tab[i][j] + (agr_tab_C[i][j] / GD.SimulationStepTime)) - sum);
//                }
//            }
//        }
//
//        for (int i = 0; i < grid1.Nn; i++)
//        {
//            double sum = 0;
//            for (int j = 0; j <= i - 1; j++)
//            {
//                sum += L[i][j] * Y[j];
//            }
//            Y[i] = agr_PC[i] - sum;
//        }
//
//        tab_wynik[grid1.Nn - 1] = Y[grid1.Nn - 1] / U[grid1.Nn - 1][grid1.Nn - 1];
//        int i = grid1.Nn - 2;
//        while (i >= 0)
//        {
//            double sum = 0;
//            for (int j = i + 1; j < grid1.Nn; j++)
//            {
//                sum += U[i][j] * tab_wynik[j];
//            }
//            tab_wynik[i] = 1 / U[i][i] * (Y[i] - sum);
//            i--;
//        }
//        cout << endl;
//        cout << "wartosc T dla: " << p << endl;
//        for (int i = 0; i < grid1.Nn; i++)
//        {
//            cout << tab_wynik[i] << "  ";
//            if (i % 4 == 3)
//            {
//
//                cout << endl;
//            }
//        }
//
//        cout << endl;
//        int zmien = 3;
//        for (int g = 0; g < grid1.Nn; g++)
//        {
//            agr_PC[g] = agr_P[g] + tab_wynik[g] * agr_C[g] / GD.SimulationStepTime;
//            cout << agr_PC[g] << "      ";
//        }
//        cout << endl;
//    }
//
//    return NULL;
//}
