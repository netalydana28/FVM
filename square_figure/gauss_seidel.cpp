#include <bits/stdc++.h>  
using namespace std;


const int center_num=3;
const double e=0.0001;

double gauss_seidel_coeff[center_num][center_num];
double gauss_seidel_b[center_num];

vector<double> P_centers(center_num);
vector<double> P_centers_copy(center_num);

double find_max(vector<double> v){
    double max = v[0];
    for(int i=0; i<v.size(); i++){
        if(v[i]>max){
            max = v[i];
        }
    }
    return max;
}

void gauss_seidel(){
    // for(int i=0; i<center_num; i++){
    //     double temp = gauss_seidel_b[i];
    //     for(int j=0; j<center_num; j++){
    //         if(i!=j){
    //             temp -= gauss_seidel_coeff[i][j]*0;
    //         } 
    //     }
    //     P_centers[i] = (1/gauss_seidel_coeff[i][i]) * 0.8;
    //     P_centers_copy[i] = P_centers[i]; 
    //     cout<<i<<"$"<<endl;

    // }
    // for(int i=0; i<center_num; i++){
    //     cout<<P_centers[i]<<" ";
    //     if(i%9 == 8){
    //         cout<<endl;
    //     }
    // }
    // while(abs(find_max(P_centers_copy) - find_max(P_centers)) > e){
    for(int m=0; m<4; m++){
    
        for(int i=0; i<center_num; i++){
            P_centers_copy[i] = (gauss_seidel_b[i] / gauss_seidel_coeff[i][i]);

            // double temp = gauss_seidel_b[i];
            for(int j=0; j<center_num; j++){
                // if(i!=j){
                //     temp -= gauss_seidel_coeff[i][j]*P_centers[i];
                // } 
                if (j == i){
                    continue;
                }

                P_centers_copy[i] = P_centers_copy[i] - ((gauss_seidel_coeff[i][j] / gauss_seidel_coeff[i][i]) * P_centers[j]);

                P_centers[i] = P_centers_copy[i];
            }
            // P_centers[i] = (1/gauss_seidel_coeff[i][i]) * temp;
            // P_centers_copy[i] = P_centers[i]; 

           
        }
        cout<<m<<"it"<<endl;
            for(int k=0; k<3; k++){
                cout<<P_centers[k]<<" ";
            }
            cout<<endl;

    }
}

void jacobi(){
    double to_check = 1;
    while(abs(to_check - find_max(P_centers))>e){
    // for(int s=0; s<4; s++){
    
        to_check = find_max(P_centers);
        for(int m=0; m<center_num; m++){
            P_centers_copy[m] = P_centers[m];
        }
        for(int i=0; i<center_num; i++){
            P_centers[i] = gauss_seidel_b[i];
            for(int j=0; j<center_num; j++){
                if (j == i){
                    continue;
                }
                P_centers[i] -= (gauss_seidel_coeff[i][j] * P_centers_copy[j]);
                cout<<P_centers_copy[j]<<" "<<j<<endl;
            }
            P_centers[i] /= gauss_seidel_coeff[i][i];

            // P_centers[i] = P_centers_copy[i];
        }
        for(int i=0; i<center_num; i++){
            cout<<P_centers[i]<<" ";
        }
        cout<<endl;
    }
}

int main(){
        //Centers
    for(int i=0; i<center_num; i++){
        P_centers.at(i)=0;
        P_centers_copy.at(i)=0;

    } 

        for (int i = 0; i < 3; i++)

    {

        for (int j = 0; j < 3; j++)

        {

            cout << "Enter values no :(" << i << ", " << j << ") ";

            cin >> gauss_seidel_coeff[i][j];

        }

    }

    cout << "\nEnter Values to the right side of equation\n";

    for (int i = 0; i < 3; i++)

    {

        cout << "Enter values no :(" << i <<") ";

        cin >> gauss_seidel_b[i];

    }

    jacobi();
    // gauss_seidel();
    for(int i=0; i<center_num; i++){
        cout<<P_centers.at(i)<<" ";
        // if(i%3 ==2){
        //     cout<<endl;
        // }

    } 
    return 0;
}