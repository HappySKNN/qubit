#include "qubit.h"

using namespace std;

Qubit::Qubit() {
    random_device r_d;
    mt19937 dist(r_d());
    uniform_real_distribution<double> random(0, 1);

    this->a_0 = random(dist);
    this->a_1 = sqrt(1 - pow(this->a_0, 2));

    Calculate_Chances();

    if (Check_Correct()) {
        cout << "Qubit was created:" << endl;
        Print_Qubit();
    } else {
        cout << "Qubit was created with errors:" << endl;
        Print_Qubit();
        exit(1);
    }
}

Qubit::Qubit(double a_0, double a_1) {
    this->a_0 = a_0;
    this->a_1 = a_1;

    Calculate_Chances();

    if (Check_Correct()) {
        cout << "Qubit was created:" << endl;
        Print_Qubit();
    } else {
        cout << "Qubit was created with errors:" << endl;
        Print_Qubit();
        exit(1);
    }
}

void Qubit::Print_Qubit() const {
    cout << this->p_0 << " -> (1, 0) = |0>" << endl
         << this->p_1 << " -> (0, 1) = |1>" << endl << endl;
}

void Qubit::Calculate_Chances() {
    this->p_0 = pow(this->a_0, 2);
    this->p_1 = pow(this->a_1, 2);
}

bool Qubit::Check_Correct() const {
    return abs(1 - (this->p_0 + this->p_1)) <= EPS;
}

void Qubit::Hadamard_Multiply() {
    vector<vector<double>> hadamard_matrix = {{1, 1},
                                              {1, -1}};

    double temp_a_0 = this->a_0, temp_a_1 = this->a_1;

    this->a_0 = (1 / sqrt(2)) * (hadamard_matrix[0][0] * temp_a_0 + hadamard_matrix[0][1] * temp_a_1);
    this->a_1 = (1 / sqrt(2)) * (hadamard_matrix[1][0] * temp_a_0 + hadamard_matrix[1][1] * temp_a_1);


    Calculate_Chances();

    if (Check_Correct()) {
        cout << "Qubit's rotated basis:" << endl;
        Print_Qubit();
    } else {
        cout << "Qubit's rotated basis was created with errors:" << endl;
        Print_Qubit();
        exit(2);
    }
}

double Qubit::Qubit_Observe() const {
    random_device r_d;
    mt19937 dist(r_d());
    uniform_real_distribution<double> random(0, 1);

    double choice = random(dist);

    if (choice <= this->p_0) {
        cout << "State: a0 = " << this->a_0 << "; Chance: " << this->p_0 << endl << endl;
        return this->a_0;
    } else {
        cout << "State: a1 = " << this->a_1 << "; Chance: " << this->p_1 << endl << endl;
        return this->a_1;
    }
}

void Qubit::Get_Matrix(const string &filename) {
    ifstream in(filename);

    vector<vector<int>> vec;

    int temp;
    for (string line; getline(in, line);) {
        istringstream iss(line);
        vector<int> vec_temp;
        while (!iss.eof()) {
            iss >> temp;
            vec_temp.push_back(temp);
        }

        vec_temp.push_back(0);
        vec_temp.push_back(vec_temp[vec_temp.size() - 2] ^ vec_temp[vec_temp.size() - 1]);
        vec.push_back(vec_temp);

        vec_temp[vec_temp.size() - 2] = 1;
        vec_temp[vec_temp.size() - 1] = vec_temp[vec_temp.size() - 3] ^ vec_temp[vec_temp.size() - 2];
        vec.push_back(vec_temp);
    }

    vector<int> dense_matrix(8);

    for (auto &i: vec) {
        stringstream ss;
        for (int j = 0; j < i.size() - 3; ++j) {
            ss << i[j];
        }
        ss << i[i.size() - 1];
        dense_matrix[&i - &vec[0]] = stoi(ss.str(), nullptr, 2);
    }

    ofstream out("result");
    for (int i = 0; i < dense_matrix.size(); ++i) {
        for (int j = 0; j < dense_matrix.size(); ++j) {
            if (dense_matrix[i] == j) {
                out << "1 ";
            } else {
                out << "0 ";
            }
        }
        out << "\n";
    }
    out.close();
}

vector<int> Qubit::Transpose_Matrix(const std::vector<int> &vec) {
    vector<int> result(8);

    for (int i = 0; i < vec.size(); ++i) {
        result[vec[i]] = i;
    }

    return result;
}

vector<int> Qubit::Multiply_Matrix(const std::vector<int> &vec) {

}