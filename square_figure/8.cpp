#include <bits/stdc++.h>  
using namespace std; 

const int n=11; 
double dx=0.25;
double P[n][n];
double P_copy[n][n];
double e = 0.001;

int main(){


    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(i==0){
                P[i][j]=1;
                P_copy[i][j]=1;

            }else{
                P[i][j]=0;
                P_copy[i][j]=0;

            }
        }
    }
    for(int k=0; k<10; k++){
        for(int i=1; i<n-1; i++){
            for(int j=1; j<n-1; j++){
                if(i % 2 == 1 and j % 2 == 1){
                    P[i][j] = (P[i-1][j-1] + P[i-1][j+1] + P[i+1][j+1] + P[i+1][j-1]) / 4;
                }else if(i % 2 == 0 and j % 2 == 0){
                    P[i][j] = (P[i-2][j] + P[i][j+2] + P[i+2][j] + P[i][j-2]) / 4;
                }else if(i % 2 == 0 and j % 2 == 1){
                    P[i][j] = (P[i][j-1] + P[i][j+1]) / 2;
                }else if(i % 2 == 1 and j % 2 == 0){
                    P[i][j] = (P[i+1][j] + P[i-1][j]) / 2;
                }
            }
        }    
    }

   
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            cout<<P[j][i]<<" ";
        }
        cout<<endl;
    }

    return 0;
}