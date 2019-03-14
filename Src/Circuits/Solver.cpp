#include <iostream>
#include <complex>
#include <cmath>
#include <string>
#include <fstream>
#include <Eigen/Eigen>
#include "Component.h"

using namespace Eigen;
using namespace std;


const double PI = 3.141592;
double w;
int VS_num = 0, CS_num = 0, R_num = 0, L_num = 0, C_num = 0, nodes_num = 0;
ifstream input;
ofstream output;
vector<Component> components;


// Read The input from external file
bool read()
{
    string s;
    cout << "Enter INPUT file name: ";
    cin >> s;
    input.open(s + ".txt");
    if (input.is_open())
    {
        input >> s;
        w = stoi(s);

        while (!input.eof())
        {
            input >> s;
            if (s == "Vs")
            {
                string n1, n2, mag, phase;
                Component VScource;
                VScource.name = s;
                input >> n1 >> n2 >> mag >> phase;
        
				VScource.node1 = stoi(n1);
                VScource.node2 = stoi(n2);
                
				double magn = stod(mag);
                double ph = stod(phase);
                
				VScource.voltage = complex<double>(magn * cos(ph * PI / 180),
                    magn * sin(ph * PI / 180));

                components.push_back(VScource);
                VS_num++;
            }
            else if (s == "Cs")
            {
                string n1, n2, mag, phase;
                Component SScource;
                SScource.name = s;
                
				input >> n1 >> n2 >> mag >> phase;
                
				SScource.node1 = stoi(n1);
                SScource.node2 = stoi(n2);
                
				double magn = stod(mag);
                double ph = stod(phase);
                
				SScource.current = complex<double>(magn * cos(ph * PI / 180),
                    magn * sin(ph * PI / 180));
                components.push_back(SScource);

                CS_num++;
            }

            else if (s[0] == 'R')
            {
                string n1, n2, mag;
                Component R;
                R.name = s;
                input >> n1 >> n2 >> mag;
                R.node1 = stoi(n1);
                R.node2 = stoi(n2);

                float magn = stof(mag);
                R.Z = complex<double>(magn, 0);
                
				R.Y = pow(R.Z, -1);
                R.factor = 0;
                components.push_back(R);
                R_num++;
            }
            else if (s[0] == 'C')
            {
                string n1, n2, factor;
                Component C;
                C.name = s;
                
				input >> n1 >> n2 >> factor;
                
				C.node1 = stoi(n1);
                C.node2 = stoi(n2);
                double magn = stod(factor);
                
				C.Z = complex<double>(0, -1 / (w * magn));
                C.Y = pow(C.Z, -1);
                
				C.factor = magn;
                components.push_back(C);
                C_num++;
            }
            else if (s[0] == 'L')
            {
                string n1, n2, factor;
                Component L;
                L.name = s;
                
				input >> n1 >> n2 >> factor;
                
				L.node1 = stoi(n1);
                L.node2 = stoi(n2);
                
				double magn = stod(factor);
                
				L.Z = complex<double>(0, w * magn);
                L.Y = pow(L.Z, -1); //(0, -1/(w*magn));
                L.factor = magn;
                
				components.push_back(L);
                L_num++;
            }
        }
    }
    else
    {
        cout << "Wrong name, Enter a valid name\n";
        return false;
    }
    return true;
}

int main()
{
    if (read())
        ;
    else
    {
        while (true)
        {
            cout << "Please Enter a valid name\n";
            if (read())
                break;
        }
    }

    for (int i = 0; i < components.size(); ++i)
    {
        if (components[i].node1 > nodes_num)
            nodes_num = components[i].node1;
        if (components[i].node2 > nodes_num)
            nodes_num = components[i].node2;
    }

    int Eqs = nodes_num;
    for (int i = 0; i < components.size(); ++i)
    {
        if (components[i].name == "Vs")
        {
            Eqs++;
            components[i].VS_num = Eqs;
        }
    }

    Eigen::MatrixXcd Ys(Eqs, Eqs);
    Eigen::MatrixXcd VS(Eqs, 1);
    Eigen::MatrixXcd IS(Eqs, 1);

    for (int i = 0; i < Eqs; ++i)
        for (int j = 0; j < Eqs; ++j)
            Ys(i, j).imag(0), Ys(i, j).real(0);

	for (int i = 0; i < Eqs; ++i)
    {
        VS(i, 0).imag(0);
        VS(i, 0).real(0);

        IS(i, 0).imag(0);
        IS(i, 0).real(0);
    }
    
	vector<complex<double> > Vnodes(nodes_num);
    for (int i = 0; i < components.size(); ++i)
    {
        // first of all, if we have a voltage source
        if (components[i].node1 == 0 && components[i].name == "Vs")
        {
            // Vnodes[components[i].node2] = components[i].voltage;
            Ys(components[i].node2 - 1, components[i].VS_num - 1) += complex<double>(-1, 0);
            Ys(components[i].VS_num - 1, components[i].node2 - 1) += complex<double>(-1, 0);
            VS(components[i].VS_num - 1, 0) += components[i].voltage;
        }
        else if (components[i].node2 == 0 && components[i].name == "Vs")
        {
            Ys(components[i].node1 - 1, components[i].VS_num - 1) += complex<double>(1, 0);
            Ys(components[i].VS_num - 1, components[i].node1 - 1) += complex<double>(1, 0);
            VS(components[i].VS_num - 1, 0) += components[i].voltage;
        }
        else if (components[i].node1 != 0 && components[i].node2 != 0
            && components[i].name == "Vs") // same as previus but we have two additional equations
        {
            Ys(components[i].node1 - 1, components[i].VS_num - 1) += complex<double>(1, 0);
            Ys(components[i].node2 - 1, components[i].VS_num - 1) += complex<double>(-1, 0);
            Ys(components[i].VS_num - 1, components[i].node1 - 1) += complex<double>(1, 0);
            Ys(components[i].VS_num - 1, components[i].node2 - 1) += complex<double>(-1, 0);
            VS(components[i].VS_num - 1, 0) += components[i].voltage;
        }

        // If we have current sources
        else if (components[i].name == "Cs" && components[i].node1 == 0)
        {
            VS(components[i].node2 - 1, 0) -= components[i].current;
        }
        else if (components[i].name == "Cs" && components[i].node2 == 0)
        {
            VS(components[i].node1 - 1, 0) += components[i].current;
        }
        else if (components[i].name == "Cs" && components[i].node1 != 0 && components[i].node2 != 0)
        {
            VS(components[i].node2 - 1, 0) -= components[i].current;
            VS(components[i].node1 - 1, 0) += components[i].current;
        }

        // Now if we have passive components
        else if (components[i].name[0] == 'R' || components[i].name[0] == 'C'
            || components[i].name[0] == 'L')
        {
            if (components[i].node1 == 0)
            {
                Ys(components[i].node2 - 1, components[i].node2 - 1) += components[i].Y;
            }
            else if (components[i].node2 == 0)
            {
                Ys(components[i].node1 - 1, components[i].node1 - 1) += components[i].Y;
            }
            else if (components[i].node1 != 0 && components[i].node2 != 0)
            {
                Ys(components[i].node1 - 1, components[i].node1 - 1) += components[i].Y;
                Ys(components[i].node2 - 1, components[i].node1 - 1) -= components[i].Y;
                Ys(components[i].node2 - 1, components[i].node2 - 1) += components[i].Y;
                Ys(components[i].node1 - 1, components[i].node2 - 1) += components[i].Y;
            }
        }
    }

    IS += Ys.jacobiSvd(ComputeThinU | ComputeThinV).solve(VS);
    cout << "\n";

    cout << Ys << "\n\n\n";
    cout << VS << "\n\n\n";
    cout << IS;
    cout << "\n";


    for (int i = 0; i < Eqs; i++)
    {

        IS(i, 0).real(int(IS(i, 0).real() * 1000) / 1000.0);
        IS(i, 0).imag(int(IS(i, 0).imag() * 1000) / 1000.0);
    }
    int i = 0;
    for (; i < nodes_num; i++)
    {
        cout << "V(" << i + 1 << ")"
             << " " << abs(IS(i, 0)) << " " << (arg(IS(i, 0)) * (180 / 3.141592654));
        cout << endl;
    }
    for (int j = 0; j < nodes_num; j++)
    {
        if (components[j].name == "Vs")
        {
            cout << "I(" << components[j].node1 << " " << components[j].node2 << ")"
                 << " " << abs(IS(i, 0)) << " " << (arg(IS(i, 0)) * (180 / 3.141592654));
            i++;
            cout << endl;
        }
        if (components[j].name == "Is")
        {
            cout << "I(" << components[j].node2 << " " << components[j].node1 << ")"
                 << " " << abs(components[j].current) << " "
                 << (arg(components[j].current) * (180 / 3.141592654));
            cout << endl;
        }
        if (components[j].name[0] == 'R' || components[j].name[0] == 'L'
            || components[j].name[0] == 'C')
        {
            if (components[j].node1 != 0 && components[j].node1 != 0)
            {
                complex<double> I
                    = ((IS(components[j].node1 - 1, 0) - IS(components[j].node2 - 1, 0))
                        / components[j].Z);
                cout << "I(" << components[j].node1 << " " << components[j].node2 << ")"
                     << " " << abs(I) << " " << (arg(I) * (180 / 3.141592654));
                cout << endl;
            }
            else
            {
                if (components[j].node1 == 0)
                {
                    complex<double> I = ((complex<double>(0, 0) - IS(components[j].node2 - 1, 0))
                        / components[j].Z);
                    cout << "I(" << components[j].node1 << " " << components[j].node2 << ")"
                         << " " << abs(I) << " " << (arg(I) * (180 / 3.141592654));
                    cout << endl;
                }
                else if (components[j].node2 == 0)
                {
                    complex<double> I = ((IS(components[j].node1 - 1, 0)) / components[j].Z);
                    cout << "I(" << components[j].node1 << " " << components[j].node2 << ")"
                         << " " << abs(I) << " " << (arg(I) * (180 / 3.141592654));
                    cout << endl;
                }
            }
        }
    }
    system("PAUSE");
    return 0;
}