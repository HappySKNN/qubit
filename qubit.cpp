#include "qubit.h"
#include <cmath>

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

void Qubit::Solve_Function(const int &a, const int &m) {
    string filename = Qubit::Get_Arithmetic_Table(a, m);
    Qubit::Get_Matrix(filename);
}

void Qubit::Solve_Function(const std::string &filename) {
    Qubit::Get_Matrix(filename);
}

void Qubit::Get_Matrix(const string &filename) {
    vector<vector<int>> vec = Get_Vector_From_Matrix(filename);

    int f_count = static_cast<int>((vec[0].size() + 1) / 2);
    int x_count = f_count;
    int lines_count = static_cast<int>(pow(pow(2, f_count), 2));
    vector<vector<int>> vec_table(lines_count, vector<int>(f_count * 4));

    if (vec[0].size() % 2) {
        lines_count /= 2;
        vec_table.resize(lines_count, vector<int>(f_count * 4 - 3));
        f_count /= 2;
    }

    for (int i = 0; i < lines_count; ++i) {
        for (int j = 0; j < x_count; ++j) {
            vec_table[i][j] = vec[i / pow(2, f_count)][j];
        }

        for (int j = 0; j < f_count; ++j) {
            vec_table[i][x_count + j] = vec[i / pow(2, f_count)][x_count + j];
        }

        int temp = i % (int) pow(2, f_count);
        for (int j = 0; j < f_count; ++j) {
            vec_table[i][x_count + f_count * 2 - j - 1] = temp & 1;
            temp >>= 1;
        }

        for (int j = 0; j < f_count; ++j) {
            vec_table[i][x_count + f_count * 3 - j - 1] =
                    vec_table[i][x_count + f_count - j - 1] ^ vec_table[i][x_count + f_count * 2 - j - 1];
        }
    }

    vector<int> dense_matrix(lines_count);

    for (auto &line: vec_table) {
        stringstream ss;

        for (int j = 0; j < x_count; ++j) {
            ss << line[j];
        }

        for (int j = x_count + f_count * 2; j < x_count + f_count * 3; ++j) {
            ss << line[j];
        }

        dense_matrix[&line - &vec_table[0]] = stoi(ss.str(), nullptr, 2);
    }

    vector<int> result_matrix = Dense_Matrices_Multiply(dense_matrix, Transpose_Matrix(dense_matrix));

    for (int i = 0; i < result_matrix.size(); ++i) {
        if (result_matrix[i] != i) {
            cout << "Matrix is not unitary" << endl;
            exit(3);
        }
    }

    ofstream out("../result_function_matrix.txt");
    for (int i = 0; i < dense_matrix.size(); ++i) {
        for (int j = 0; j < dense_matrix.size(); ++j) {
            if (j < dense_matrix.size() - 1) {
                if (dense_matrix[i] == j) {
                    out << "1 ";
                } else {
                    out << "0 ";
                }
            } else {
                if (dense_matrix[i] == j) {
                    out << "1";
                } else {
                    out << "0";
                }
            }
        }
        if (i < dense_matrix.size() - 1) {
            out << endl;
        }
    }
    out.close();

    cout << "Matrix was created" << endl;
}

vector<int> Qubit::Transpose_Matrix(const std::vector<int> &vec) {
    vector<int> result(vec.size());

    for (int i = 0; i < vec.size(); ++i) {
        result[vec[i]] = i;
    }

    return result;
}

vector<int> Qubit::Dense_Matrices_Multiply(const std::vector<int> &matrix_1, const std::vector<int> &matrix_2) {

    vector<int> result(matrix_1.size());

    for (int i = 0; i < matrix_1.size(); ++i) {
        result[i] = matrix_2[matrix_1[i]];
    }

    return result;
}

vector<vector<int>> Qubit::Multiply_Matrices(vector<vector<int>> matrix_1, vector<vector<int>> matrix_2) {
    vector<vector<int>> result(matrix_1.size(), vector<int>(matrix_2[0].size()));

    for (int i = 0; i < matrix_1.size(); ++i) {
        for (int j = 0; j < matrix_2[0].size(); ++j) {
            for (int k = 0; k < matrix_2.size(); ++k) {
                result[i][j] += matrix_1[i][k] * matrix_2[k][j];
            }
        }
    }

    return result;
}

string Qubit::Get_Arithmetic_Table(const int &a, const int &m) {
    int f_count = ceil(log2(m));
    int lines_count = static_cast<int>(pow(2, f_count));

    string filename = "../arithmetic_table.txt";
    ofstream out(filename);

    for (int i = 0; i < lines_count; ++i) {
        vector<int> to_write(f_count * 2);

        int temp = i;
        for (int j = 0; j < f_count; ++j) {
            to_write[f_count - j - 1] = temp & 1;
            temp >>= 1;
        }

        temp = ((int) a * i % m);
        for (int j = 0; j < f_count; ++j) {
            to_write[f_count * 2 - j - 1] = temp & 1;
            temp >>= 1;
        }

        for (int j = 0; j < f_count * 2; ++j) {
            if (j < f_count * 2 - 1) {
                out << to_write[j] << " ";
            } else {
                out << to_write[j];
            }

        }


        if (i < lines_count - 1) {
            out << "\n";
        }
    }

    out.close();

    return filename;
}

void Qubit::Get_Qubit(const string &filename) {
    Solve_Function(filename);
    vector<vector<int>> const matrix_vector = Get_Vector_From_Matrix("../result_function_matrix.txt");

    vector<vector<int>> boolean_matrix_vector = Get_Vector_From_Matrix(filename);
    for (auto &line: boolean_matrix_vector) {
        line.back() = 0;
    }

    vector<int> linear_vector(matrix_vector.size(), 0);

    for (auto &line: boolean_matrix_vector) {
        int temp = 0;
        for (int i = 0; i < line.size(); ++i) {
            temp += line[i] * (int) pow(2, line.size() - i - 1);
        }
        linear_vector[temp] = 1;
    }

    vector<vector<int>> matrix_from_linear_vector(linear_vector.size(), vector<int>(linear_vector.size()));
    matrix_from_linear_vector[0] = linear_vector;

    vector<vector<int>> const result_vector = Multiply_Matrices(matrix_from_linear_vector, matrix_vector);

    cout << "State: 1/2 * (";
    for (int i = 0; i < result_vector.size(); ++i) {
        if (result_vector[0][i] == 1) {
            cout << "|" << Convert_Decimal_To_Binary(i, int(ceil(log(result_vector.size())))) << ">";
            if (i != result_vector.size() - 1) {
                cout << " + ";
            }
        }
    }
    cout << ")" << endl;
}

void Qubit::Get_Qubit(const vector<int> &boolean_vector) {

}

std::vector<std::vector<int>> Qubit::Get_Vector_From_Matrix(const string &filename) {
    ifstream in(filename);
    vector<vector<int>> vec;

    for (string line; getline(in, line);) {
        istringstream iss(line);
        int temp;
        vector<int> vec_temp;

        while (!iss.eof()) {
            iss >> temp;
            vec_temp.push_back(temp);
        }

        vec.push_back(vec_temp);
    }

    return vec;
}

string Qubit::Convert_Decimal_To_Binary(const int &number, const int &length) {
    string result;

    int temp = number;
    while (temp > 0) {
        result += to_string(temp & 1);
        temp >>= 1;
    }

    while (result.size() < length) {
        result += "0";
    }

    for (int i = 0; i < result.size() / 2; ++i) {
        swap(result[i], result[result.size() - i - 1]);
    }

    return result;
}

std::vector<std::vector<int>>
Qubit::Tensor_Multiply(const vector<std::vector<int>> &matrix_1, const vector<std::vector<int>> &matrix_2) {
    vector<vector<int>> result(matrix_1.size() * matrix_2.size(), vector<int>(matrix_1[0].size() * matrix_2[0].size()));

    for (int i = 0; i < matrix_1.size(); ++i) {
        for (int j = 0; j < matrix_1[0].size(); ++j) {
            for (int k = 0; k < matrix_2.size(); ++k) {
                for (int l = 0; l < matrix_2[0].size(); ++l) {
                    result[i * matrix_2.size() + k][j * matrix_2[0].size() + l] = matrix_1[i][j] * matrix_2[k][l];
                }
            }
        }
    }

    return result;
}

void Qubit::Deutsch_Jozsa_Algorithm(const string &filename) {
    Solve_Function(filename);
    vector<vector<int>> const matrix_vector = Get_Vector_From_Matrix("../result_function_matrix.txt");

    //vector<vector<int>> boolean_matrix_vector = Get_Vector_From_Matrix(filename);

    vector<vector<int>> Hadamard_Matrix = {{1, 1},
                                           {1, -1}};

    vector<vector<int>> multiply_hadamard_matrix = Tensor_Multiply(Hadamard_Matrix, Hadamard_Matrix);

    multiply_hadamard_matrix = Tensor_Multiply(multiply_hadamard_matrix, Hadamard_Matrix);

    vector<vector<int>> vec(matrix_vector.size(), vector<int>(matrix_vector.size(), 1));
    for (int i = 0; i < vec.size(); ++i) {
        if (i % 2)
            vec[0][i] *= -1;
    }

    vector<vector<int>> result_vector = Multiply_Matrices(vec, matrix_vector);

    result_vector = Multiply_Matrices(result_vector, multiply_hadamard_matrix);

    cout << "State: (";
    for (int i = 0; i < result_vector.size(); ++i) {
        if (result_vector[0][i] != 0) {
            cout << "|" << Convert_Decimal_To_Binary(i, int(ceil(log(result_vector.size())))) << ">";
            if (i != result_vector.size() - 1) {
                cout << " + ";
            }
        }
    }
    cout << ")" << endl;
}
