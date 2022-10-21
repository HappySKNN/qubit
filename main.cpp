#include "qubit.h"

using namespace std;

void Task_1() {
    Qubit a;
    Qubit b;

    a.Qubit_Observe();
    b.Hadamard_Multiply();
    b.Qubit_Observe();
}

void Task_2() {
    random_device r_d;
    mt19937 dist(r_d());
    uniform_int_distribution<int> random(0, 3);

    vector<Qubit> vec = {Qubit(0, 1), Qubit(1, 0), Qubit(0, 1), Qubit(1, 0)};
    vec[2].Hadamard_Multiply();
    vec[3].Hadamard_Multiply();

    string result;

    for (int i = 0; i < 8;) {
        int j = random(dist);
        double res;

        if (!random(dist)) {
            res = vec[j].Qubit_Observe();
        } else {
            vec[j].Hadamard_Multiply();
            res = vec[j].Qubit_Observe();
            vec[j].Hadamard_Multiply();
        }

        if (1 - res <= EPS || res <= EPS) {
            cout << "OK" << endl;
            if (j == 0 || j == 2) {
                result += "0";
            } else {
                result += "1";
            }
            i++;
        } else {
            cout << "REPEAT" << endl;
        }
    }
    cout << result;
}

int main() {
    //Task_1();
    //Task_2();
    //Qubit::Get_Matrix("../initial_table.txt");
    Qubit::Calculate_Function(3, 5);
    return 0;
}
