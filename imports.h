#include <iostream>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>
using namespace std;

    //variables & constants//
    int i, j, N, n, num, k, l, num_a;
    float acid_table[20][20];
    char bond_table[20][20];
    char acids[20];
    const double b0 = 1.32;
    const double kb = 8.617333262135;
    string line, seq_v0, probability, letter;
    int key = 'X';
    const float c_0 = 1.33;
    const float c_1 = 1.44;
    const float d_min = 1.55;

    map<char, int> acid_num;
    map<char, int> atoms_per_amino;

    //create a vector of file names//
    typedef vector<string> stringvec;

int get_maps_and_tables() {
    //hydrogen bond probability chart//
	ifstream prob("AminoAcidTable0.txt");
	for (i = 0; i < 20; ++i) {
        for (j = 0; j < 20; ++j) {
            prob >> acid_table[i][j];
            //cout << acid_table[i][j] << " \n"[j == 20-1];
        }
	}
    //hydrogen bond type chart//
    ifstream bond_chart("sidechains.txt");
	for (i = 0; i < 20; ++i) {
        for (j = 0; j < 20; ++j) {
            bond_chart >> bond_table[i][j];
            //cout << bond_table[i][j] << " \n"[j == 20-1];
        }
	}

    //amino acid indices//
    acid_num['A']=0;
    acid_num['R']=1;
    acid_num['N']=2;
    acid_num['D']=3;
    acid_num['C']=4;
    acid_num['Q']=5;
    acid_num['E']=6;
    acid_num['G']=7;
    acid_num['H']=8;
    acid_num['I']=9;
    acid_num['L']=10;
    acid_num['K']=11;
    acid_num['M']=12;
    acid_num['F']=13;
    acid_num['P']=14;
    acid_num['S']=15;
    acid_num['T']=16;
    acid_num['W']=17;
    acid_num['Y']=18;
    acid_num['V']=19;

    //number of atoms per amino acid index//
    atoms_per_amino[0]=10;
    atoms_per_amino[1]=11;
    atoms_per_amino[2]=8;
    atoms_per_amino[3]=8;
    atoms_per_amino[4]=6;
    atoms_per_amino[5]=18;
    atoms_per_amino[6]=9;
    atoms_per_amino[7]=4;
    atoms_per_amino[8]=10;
    atoms_per_amino[9]=8;
    atoms_per_amino[10]=16;
    atoms_per_amino[11]=9;
    atoms_per_amino[12]=8;
    atoms_per_amino[13]=11;
    atoms_per_amino[14]=7;
    atoms_per_amino[15]=6;
    atoms_per_amino[16]=7;
    atoms_per_amino[17]=14;
    atoms_per_amino[18]=12;
    atoms_per_amino[19]=14;
    return 0;
}

//returns the end to end distance of the polymer chain//
double R (vector<int> v1, vector<int> v2, vector<int> v3) {
    double r;
    r = sqrt(pow(v1.begin() - v1.end(), 2) + pow(v2.begin() - v2.end(), 2) + pow(v3.begin() - v3.end(), 2) );
    return r;
}

double Entropy (double var1, double var2, double var3) {
    double var_total;
    var_total = var1 + var2 + var3;
    return var_total;
}

void read_directory(const string& name, stringvec& v)
{
    DIR* dirp = opendir(name.c_str());
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL) {
        v.push_back(dp->d_name);
    }
    closedir(dirp);
}

int f_elec_hi(vector<vector <float> > var4) {
    return 0;
}

int F_desolv(vector<vector <float> > var5) {
    return 0;
}

int F_desolv_low(vector<vector <float> > var5) {
    return 0;
}

int F_desolv_high(vector<vector <float> > var6) {
    return 0;
}

int Ef_poly(vector<vector <float> > var7) {
    return 0;
}

int E_elec(vector<vector <float> > var8) {
    return 0;
}

int Efa_elec(vector<vector <float> > var9) {
    return 0;
}