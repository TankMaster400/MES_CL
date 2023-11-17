/*
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
    cout << endl;

    for (int j = 0; j < 4; j++)
    {
        for (int g = 0; g < 4; g++)
        {
            cout << "   " << H.H_k[j][g];
        }
        cout << endl;
    }

    */

/*
struct matrixH
{
    double** dNdx;
    double** dNdy;
    double detJ[4];
    double tk[4][2][2];
    double H[4][4][4];
    double H_k[4][4];
    matrixH(double ID_xy[2][4], struct element_uni el)
    {
        for (int i = 0; i < el.N; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                H_k[i][j] = 0;
            }
        }

        for (int i = 0; i < 4; i++)
        {
            tk[i][0][0] = el.dNdKsi[i][0] * ID_xy[0][0] + el.dNdKsi[i][1] * ID_xy[0][1] + el.dNdKsi[i][2] * ID_xy[0][2] + el.dNdKsi[i][3] * ID_xy[0][3];
            tk[i][0][1] = el.dNdKsi[i][0] * ID_xy[1][0] + el.dNdKsi[i][1] * ID_xy[1][1] + el.dNdKsi[i][2] * ID_xy[1][2] + el.dNdKsi[i][3] * ID_xy[1][3];
            tk[i][1][0] = el.dNdEta[i][0] * ID_xy[0][0] + el.dNdEta[i][1] * ID_xy[0][1] + el.dNdEta[i][2] * ID_xy[0][2] + el.dNdEta[i][3] * ID_xy[0][3];
            tk[i][1][1] = el.dNdEta[i][0] * ID_xy[1][0] + el.dNdEta[i][1] * ID_xy[1][1] + el.dNdEta[i][2] * ID_xy[1][2] + el.dNdEta[i][3] * ID_xy[1][3];

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
                dNdx[i][j] = tk[i][0][0] * detJ[i] * el.dNdKsi[i][j] + tk[i][1][0] * detJ[i] * el.dNdEta[i][j];
            }
        }

        for (int i = 0; i < el.N; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                dNdy[i][j] = tk[i][0][1] * detJ[i] * el.dNdKsi[i][j] + tk[i][1][1] * detJ[i] * el.dNdEta[i][j];
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

        for (int j = 0; j < 4; j++)
        {
            for (int g = 0; g < 4; g++)
            {
                for (int i = 0; i < 4; i++)
                {
                    H_k[j][g] += H[i][j][g] * 1 * 1;
                }

            }
        }
    }
};
*/