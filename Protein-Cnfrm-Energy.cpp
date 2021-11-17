#include <iostream>
#include <map>
#include <fstream>
#include <numeric>
#include "imports.h"
#include <typeinfo>
using namespace std;

int main() {

    ifstream database;
	database.open ("pdb_seqres.txt");
    get_maps_and_tables();
    n = 0;

    //navigate to the mmCIF folder//
    const char *mmCIF_path = "/Users/jamesjennings/Desktop/projects/mmcif/mmCIF/";
    chdir(mmCIF_path);

    stringvec folder_list;
    stringvec cif_file_list;
    read_directory(".", folder_list);
    copy(folder_list.begin(), folder_list.end(), ostream_iterator<string>(cout, "\n"));

    for (i = 0; i < folder_list.size(); i++) {
        const char * cif_path = folder_list[i].c_str();
        chdir(cif_path);
        read_directory(".", cif_file_list);
        copy(cif_file_list.begin(), cif_file_list.end(), ostream_iterator<string>(cout, "\n"));
        chdir(mmCIF_path);
    }

    //conformational entropy//
    while (!database.eof() and n != 100) {
        getline(database,seq_v0);
        if (!(seq_v0[0] == '>' || count(seq_v0.begin(), seq_v0.end(), key) != 0)) {

            //backbone-backbone Matrtices//
            vector <vector <int> > Hbb (seq_v0.size(), vector<int>()); //bb energy//
            vector<float> Sbb (seq_v0.size()); //bb entropy//

            //sidechaine-sidechain Matrices//
            vector <vector <float> > M (seq_v0.size(), vector<float>()); //h-bond probability//
            vector <vector <float> > H (seq_v0.size(), vector<float>()); //sc energy * M//
            vector <vector <float> > Hsc (seq_v0.size(), vector<float>()); //sc energy//
            vector <vector <char> > H_char (seq_v0.size(), vector<char>()); //type of H-bond or absense//
            vector <float> Ssc (seq_v0.size()); //sc entropy//
            vector <int> seq_indices; //amino acid indices//

            //Van der Waals amino Matrices//
            vector <vector <float> > Evdw (seq_v0.size(), vector<float>());
            vector <vector <float> > Hvdw_kl (seq_v0.size(), vector<float>()); 
            vector <vector <float> > Wconn_kl (seq_v0.size(), vector<float>());

            //electrostatic amino matrices//
            vector <vector <float> > Hfa_elec (seq_v0.size(), vector<float>());
            vector <vector <float> > Helec_kl (seq_v0.size(), vector<float>());  

            //solvent preliminaries//
            vector <vector <float> > Gdesolv_kl (seq_v0.size(), vector<float>());
            vector <vector <float> > Fdesolv (seq_v0.size(), vector<float>());

            //solvent amino matrices//
            vector <vector <float> > Efa_solv (seq_v0.size(), vector<float>());
            vector <vector <float> > H_lk_ball_kl (seq_v0.size(), vector<float>());
            vector <vector <float> > E_lk_ball_iso_kl (seq_v0.size(), vector<float>());
            vector <vector <float> > Hdesolv_kl (seq_v0.size(), vector<float>());
            vector <vector <float> > Elk_ball_wtd (seq_v0.size(), vector<float>());
            vector <vector <float> > hdesolv_kl (seq_v0.size(), vector<float>());

            //Van der waals atom matrices//
            //vector <vector<float> > Eatr_mn (20, vector<float>());
            vector <vector <float> > Erep_mn;
            vector <vector <float> > Wconn_mn;

            vector <vector <float> > Hvdw_mn;

            //Van der Waals attractive and repulsive energy matrices//
            vector <vector <float> > E;
            vector <vector <float> > a;
            vector <vector <float> > b;
            vector <vector <float> > d;
            vector <vector <float> > m;
            
            //electrostatic atom matrices//
            vector <vector <float> > Hfa_elec_mn;
            vector <vector <float> > Helec_mn;

            //Shared terms for the solvation energy//
            vector <vector <float> > Gdesolv_mn;
            vector <vector <float> > Fdesolv_mn;

            vector <vector <float> > F_lkfrac_mn;
            vector <vector <float> > E_lk_ball_mn;
            vector <vector <float> > E_lk_ball_iso_mn;
            vector <vector <float> > Hdesolv_mn;

            vector <int> NumAk;
            vector <int> NumAl;

            //amino acid coordinates for the first amino acid in the amino acid pair//
            vector <int> x1;
            vector <int> y1;
            vector <int> z1;

            //amino acid coordinates for the second amino acid in the amino acid pair//
            vector <int> x2; 
            vector <int> y2;
            vector <int> z2;

            //polymer is in a straight line//
            for (int m = 0; m < seq_v0.size(); m++) {
                x1.push_back(m);
                y1.push_back(m);
                z1.push_back(m);

                x2.push_back(m);
                y2.push_back(m);
                z2.push_back(m);
            }

            //return the end to end distance// 
            //of the given polymer chain//
            R(x1, y1, z1);

            for (int k = 0; k < seq_v0.size(); k++) {
            seq_indices.push_back(acid_num[seq_v0[k]]);
            cout << seq_indices[k] << "\n";

            for (int l = 0; l < k; l++) {
            M[k].push_back(acid_table[seq_indices[k]][seq_indices[l]]);
            H_char[k].push_back(bond_table[seq_indices[k]][seq_indices[l]]);

                //average energies for the O--H and N--H interactions//
                if (k == l + 1 || l == k + 1) {
                    H[k].push_back(0);

                    } else {
                        if (H_char[k][l] == 'O') {
                        H[k].push_back(-0.001309 * M[k][l]);
                        } else if (H_char[k][l] == 'N') {
                        H[k].push_back(-0.0010237 * M[k][l]);
                        } else {
                        H[k].push_back(0);
                        }

                    }

                //side_chain entropy//
                if (x1[l] - x2[k] == 1 && y1[l] - y2[k] == 1 && z1[l] - z2[k] == 1) {
                Hsc[k].push_back(H[k][l]);

                } else {
                Hsc[k].push_back(H[k][l]);

                }

                //back_bone entropy//
                if (x1[l] - x2[k] == 1 && y1[l] - y2[k] == 1 && z1[l] - z2[k] == 1
                && x1[l + 1] - x2[k + 1] == 1 && y1[l + 1] - y2[k + 1] == 1 && z1[l + 1] - z2[k + 1] == 1
                && x1[l + 2] - x2[k + 2] == 1 && y1[l + 2] - y2[k + 2] == 1 && z1[l + 2] - z2[k + 2] == 1) {
                Hbb[k].push_back(1);

                } else if (x1[l] - x2[k] == 1 && y1[l] - y2[k] == 1 && z1[l] - z2[k] == 1
                && x1[l + 1] - x2[k - 1] == 1 && y1[l + 1] - y2[k - 1] == 1 && z1[l + 1] - z2[k - 1] == 1
                && x1[l + 2] - x2[k - 2] == 1 && y1[l + 2] - y2[k - 2] == 1 && z1[l + 2] - z2[k - 2] == 1) {
                Hbb[k].push_back(1);

                } else {
                Hbb[k].push_back(0);

                }

                //FIX ME I AM REEKING WITH BUGS//
                for (int i = 0; i < 20; i++) {
                for (int j = 0; j < i; j++) {

                    //attractive van der waals energy//
                    if (d[i][j] <= a[i][j]) {
                    Eatr_mn[i].push_back(- E[i][j]);

                    } else if (a[i][j] < d[i][j] <= 4.5) {
                    Eatr_mn[i].push_back(E[i][j] * (pow(a[i][j], 12.0) * pow(d[i][j], -12.0) - 2 * pow(a[i][j], 6.0) * pow(d[i][j], -6.0)));

                    } else if (4.5 <= d[i][j] <= 6.0) {
                    Eatr_mn[i].push_back(Ef_poly(d[i][j]));

                    } else if (6.0 <= d[i][j]) {
                    Eatr_mn[i].push_back(0);

                    }

                    //repulsive van der waals energy//
                    if (d[i][j] <= 0.6 * a[i][j] ) {
                    Erep_mn[i].push_back(m[i][j] * d[i][j] + b[i][j]);

                    } else if (0.6 * a[i][j] < d[i][j] <= 4.5) {
                    Eatr_mn[i].push_back(E[i][j] * (pow(a[i][j], 12.0) * pow(d[i][j], -12.0) - 2 * pow(a[i][j], 6.0) * pow(d[i][j], -6.0) + 1));

                    } else if (a[i][j] < d[i][j]) {
                    Eatr_mn[i].push_back(0);

                    }

                    //electrostatic energy//
                    if (d[i][j] < 1.45) {
                    Hfa_elec_mn[i].push_back(Efa_elec(d_min));

                    } else if (1.45 <= d[i][j] < 1.85) {
                    Hfa_elec_mn[i].push_back(f_elec_low(d[i][j]));

                    } else if (1.85 <= d[i][j] < 4.5) {
                    Hfa_elec_mn[i].push_back(E_elec(d[i][j]));

                    } else if (4.5 <= d[i][j] < 5.5) {
                    Hfa_elec_mn[i].push_back(f_elec_hi(d[i][j]));

                    } else if (5.5 <= d[i][j]) {
                    Hfa_elec_mn[i].push_back(0);

                    }

                    //solvent energy//
                    if (d[i][j] <= a[i][j] - c_0) {
                    Gdesolv_mn[i].push_back(F_desolv(a[i][j]));

                    } else if (a[i][j] - c_0 < d[i][j] <= a[i][j] + c_1) {
                    Gdesolv_mn[i].push_back(F_desolv_low(i, j, d[i][j]));
 
                    } else if (a[i][j] - c_1 < d[i][j] <= 4.5) {
                    Gdesolv_mn[i].push_back(F_desolv(d[i][j]));

                    } else if (4.5 < d[i][j] <= 6.0) {
                    Gdesolv_mn[i].push_back(F_desolv_high(d[i][j]));

                    } else if (6.0 < d[i][j]) {
                    Gdesolv_mn[i].push_back(0);
                        
                    }

                    Hdesolv_mn[i].push_back(E_lk_ball_iso_mn[i][j] + E_lk_ball_mn[i][j]);
                    Hvdw_mn[i].push_back(Eatr_mn[i][j] * Wconn_mn[i][j] + Erep_mn[i][j] * Wconn_mn[i][j]);

                }
                }
                    //not sure about this one//
                    //total isotropic solvation energy//
                    Efa_solv[l].push_back(Wconn_kl[k][l] * Gdesolv_kl[k][l]);

                    //not sure about this one either//
                    //total anisotropic solvation energy//
                    Elk_ball_wtd[l].push_back(Wconn_kl[k][l] * Hdesolv_kl[k][l]);
                    Evdw[l].push_back(Evdw_kl[k][l] * Wconn_kl[k][l]);

                    Sbb[k] = accumulate(Hbb[k].begin(),Hbb[k].end(),0);
                    Ssc[k] = accumulate(Hsc[k].begin(),Hsc[k].end(),0);
            }
            }
            Entropy(R(x1,y1,z1), accumulate(Sbb.begin(),Sbb.end(),0), accumulate(Ssc.begin(),Ssc.end(),0));
        }
        
        n++;
        if (n % 100 == 0) cout << " | " << n; // doesn't work for some reason, don't panic
    database.close();
    }
    return 0;
}
