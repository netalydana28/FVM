#include <bits/stdc++.h>  
using namespace std; 

const int n=10; 
double dx=0.5;
double P[n][n];
double P_copy[n][n];
double e = 0.0001;

vector<double> row_return(vector<vector<double> > v, int n){
    vector<double> row; 
    for(int i; i<v[n].size(); i++){
        row.push_back(v[n][i]);
    }
    return row;
}

vector<double> normal(vector<double> l1, vector<double> l2){ 
    vector <double> ans; 
    // vector <double> c; 
    // vector <double> d; 

   
    // cout << l1[0] << " " << l1[1] << " " << l2[0] << " " << l2[1] << endl;
    ans.push_back(l1[1] - l2[1]); 
    ans.push_back(l2[0] - l1[0]); 
    // d.push_back(l2[1] - l2[1]); 
    // d.push_back(l1[0] - l2[0]);  
    // ans.push_back(d[0] - c[0]); 
    // ans.push_back(d[1] - c[1]);   

    double len = pow(pow(ans[0], 2) + pow(ans[1], 2), 0.5); 
    for(int i=0; i<ans.size(); i++){
        ans[i] /= len;
    }
    
    return ans; 
}



// void fvm(vector<double> ind){ 
//     double P1; 
//     vector<double> s;
//     vector<vector<double> > points; 
//     vector<vector<double> > center, nor;
    
//     vector<double> temp;
//     temp.push_back(ind[0]);
//     temp.push_back(ind[1]-1);
//     points.push_back(temp);
//     temp.clear();
//     temp.push_back(ind[0]-1);
//     temp.push_back(ind[1]);
//     points.push_back(temp);
//     temp.clear();
//     temp.push_back(ind[0]-1);
//     temp.push_back(ind[1]+1);
//     points.push_back(temp);
//     temp.clear();
//     temp.push_back(ind[0]);
//     temp.push_back(ind[1]+1);
//     points.push_back(temp);
//     temp.clear();
//     temp.push_back(ind[0]+1);
//     temp.push_back(ind[1]);
//     points.push_back(temp);
//     temp.clear();
//     temp.push_back(ind[0]+1);
//     temp.push_back(ind[1]-1);
//     points.push_back(temp);
//     temp.clear();
    
//     for(int i=0; i<points.size(); i++){
//         cout<<points[i][0]<< " "<<points[i][1]<<endl;
//         if(i == points.size()-1){
//             //cout<<normal(row_return(points, i), row_return(points, 0));

//         }else{
//             //cout<<normal(row_return(points, i), row_return(points, i+1));
//         }
//         vector<double> temp;
//         temp = row_return(points, i);
//         cout<<temp<<endl;
//     }

// }

int main(){
    vector<double> ind;
    ind.push_back(1);
    ind.push_back(1);
    // fvm(ind);
    vector<vector<double> > center, nor; 
    vector<double> temp;
    temp.push_back(ind[0]);
    temp.push_back(ind[1]-1);
    center.push_back(temp);
    temp.clear();
    temp.push_back(ind[0]-1);
    temp.push_back(ind[1]);
    center.push_back(temp);
    temp.clear();
    temp.push_back(ind[0]-1);
    temp.push_back(ind[1]+1);
    center.push_back(temp);
    temp.clear();
    temp.push_back(ind[0]);
    temp.push_back(ind[1]+1);
    center.push_back(temp);
    temp.clear();
    temp.push_back(ind[0]+1);
    temp.push_back(ind[1]);
    center.push_back(temp);
    temp.clear();
    temp.push_back(ind[0]+1);
    temp.push_back(ind[1]-1);
    center.push_back(temp);
    temp.clear();
    for(int i=0; i<center.size(); i++){
        vector<double> temp;
        if(i == center.size()-1){
            temp = normal(row_return(center, i), row_return(center, 0));
        }
        else{
            temp = normal(row_return(center, i), row_return(center, i+1));
        }
    
        nor.push_back(temp);
        cout << "***";
        cout<<nor[i][0]<<" "<<nor[i][1]<<endl;
        temp.clear();
        cout << "***";
    }
    // }
    return 0;
}