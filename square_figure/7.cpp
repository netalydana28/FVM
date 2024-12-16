//Cell-centered fvm 
#include <bits/stdc++.h>  
using namespace std; 

const int n=19; 
double dx=0.25;
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

double find_max(double matrix1[n][n]){ 
    double max = 0; 
    for(int i=1; i<n-1; i++){ 
        for(int j=1; j<n-1; j++){ 
            if (matrix1[i][j] > max){ 
                max = matrix1[i][j]; 
            } 
        } 
    } 
    return max; 
}

vector<double> normal(vector<double> l1, vector<double> l2){ 
    vector <double> ans; 
   
    cout << l1[0] << " " << l1[1] << " " << l2[0] << " " << l2[1] << endl;
    ans.push_back(l1[1] - l2[1]); 
    ans.push_back(l2[0] - l1[0]);  

    double len = pow(pow(ans[0], 2) + pow(ans[1], 2), 0.5); 
    for(int i=0; i<ans.size(); i++){
        ans[i] /= len;
    }
    return ans; 
}

void fvm(vector<int> ind){ 
    double P1; 
    vector<double> s;
    vector<vector<int> > points; 
    vector<vector<double> > center, nor;
    
    vector<int> temp;
    temp.push_back(ind[0]-1);
    temp.push_back(ind[1]-1);
    points.push_back(temp);
    temp.clear();
    temp.push_back(ind[0]-1);
    temp.push_back(ind[1]+1);
    points.push_back(temp);
    temp.clear();
    temp.push_back(ind[0]+1);
    temp.push_back(ind[1]+1);
    points.push_back(temp);
    temp.clear();
    temp.push_back(ind[0]+1);
    temp.push_back(ind[1]-1);
    points.push_back(temp);
    temp.clear();
    // temp.push_back(ind[0]);
    // temp.push_back(ind[1]-1);
    // points.push_back(temp);
    // temp.clear();
    // temp.push_back(ind[0]-1);
    // temp.push_back(ind[1]);
    // points.push_back(temp);
    // temp.clear();
    // temp.push_back(ind[0]);
    // temp.push_back(ind[1]+1);
    // points.push_back(temp);
    // temp.clear();
    // temp.push_back(ind[0]+1);
    // temp.push_back(ind[1]);
    // points.push_back(temp);
    // temp.clear();
    // vector<int> temp;
    // temp.push_back(ind[0]-1);
    // temp.push_back(ind[1]-1);
    // points.push_back(temp);
    // temp.clear();
    // temp.push_back(ind[0]-1);
    // temp.push_back(ind[1]+1);
    // points.push_back(temp);
    // temp.clear();
    // temp.push_back(ind[0]+1);
    // temp.push_back(ind[1]+1);
    // points.push_back(temp);
    // temp.clear();
    // temp.push_back(ind[0]+1);
    // temp.push_back(ind[1]-1);
    // points.push_back(temp);
    // temp.clear();

    // if(ind[1]==1){
    //     temp.push_back(ind[0]);
    //     temp.push_back(ind[1]-1);
    // }else{
    //     temp.push_back(ind[0]);
    //     temp.push_back(ind[1]-2);
    // }  
    // points.push_back(temp);
    // s.push_back((ind[0]+1) - (ind[0]-1));
    // temp.clear();

    // if(ind[0]==1){
    //     temp.push_back(ind[0]-1);
    //     temp.push_back(ind[1]);
    // }else{
    //     temp.push_back(ind[0]-2);
    //     temp.push_back(ind[1]);
    // }
    // points.push_back(temp);
    // s.push_back((ind[1]+1) - (ind[1]-1));
    // temp.clear();
    // if(ind[1]==(n-2)){
    //     temp.push_back(ind[0]);
    //     temp.push_back(ind[1]+1);
    // }else{
    //     temp.push_back(ind[0]);
    //     temp.push_back(ind[1]+2);
    // }
    // points.push_back(temp);
    // s.push_back((ind[0]+1) - (ind[0]-1));
    // temp.clear();
    // if(ind[0]==(n-2)){
    //     temp.push_back(ind[0]+1);
    //     temp.push_back(ind[1]);
    // }else{
    //     temp.push_back(ind[0]+2);
    //     temp.push_back(ind[1]);
    // }
    // s.push_back((ind[1]+1) - (ind[1]-1));
    // points.push_back(temp);
    // temp.clear();

    for(int i=0; i<points.size(); i++){
        cout<<points[i][0]<<" "<<points[i][1]<<endl;
        cout<<P[points[i][0]][points[i][1]]<<endl;
    }


    for(int i=0; i<points.size(); i++){
        
        vector<double> temp;
        if(i == 0){
            temp.push_back(dx*(points[0][0] +points[points.size()-1][0])/2);
            temp.push_back(dx*(points[0][1] +points[points.size()-1][1])/2);
        }else{
            temp.push_back(dx*(points[i][0] +points[i-1][0])/2);
            temp.push_back(dx*(points[i][1] +points[i-1][1])/2);
        }

        center.push_back(temp);
        temp.clear();
        cout<<"Center"<<endl;
        cout<<center[i][0]<<" "<<center[i][1]<<endl;
    }

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
    //length 
    for(int i=0; i<center.size(); i++){
        if(i == center.size()-1){
            s.push_back(pow((pow((center[i][0] - center[0][0]), 2) + pow((center[i][1] - center[0][1]), 2)) , 0.5));
        }else{
            s.push_back(pow((pow((center[i][0] - center[i+1][0]), 2) + pow((center[i][1] - center[i+1][1]), 2)) , 0.5));
        }
        cout<<i<< " "<<s[i]<<endl;
    } 
     
    double numerator = 0; 
    double denominator = 0; 
    for(int i=0; i<points.size(); i++){
        cout<<"x: "<<points[i][0]<< " y: "<<points[i][1]<<endl;
        if(((points[i][0] - ind[0]) == 0) and ((points[i][1] - ind[1]) == 0)){
            continue;
        }else if((points[i][0] - ind[0]) == 0){
            numerator += (P[points[i][0]][points[i][1]]*nor[i][1]*s[i])/((points[i][1] - ind[1])*dx);
            cout<<"1: "<<P[points[i][0]][points[i][1]]<<" "<<nor[i][1]<<" "<<s[i]<<" "<<(points[i][1] - ind[1])*dx<<endl;
            denominator += (s[i]*nor[i][1])/((points[i][1] - ind[1])*dx);
        }else if((points[i][1] - ind[1]) == 0){
            cout<<"2: "<<P[points[i][0]][points[i][1]]<<" "<<nor[i][0]<<" "<<s[i]<<" "<<(points[i][0] - ind[0])*dx<<endl;
            numerator += (P[points[i][0]][points[i][1]]*nor[i][0]*s[i])/((points[i][0] - ind[0])*dx);
            denominator += (s[i]*nor[i][0])/((points[i][0] - ind[0])*dx);
        }else{
            cout<<"3: "<<P[points[i][0]][points[i][1]]<<" "<<nor[i][0]<<" "<<s[i]<<" "<<(points[i][0] - ind[0])*dx<<endl;
            cout<<"4: "<<P[points[i][0]][points[i][1]]<<" "<<nor[i][1]<<" "<<points[i][1]<<" "<<ind[1]*dx<<endl;
            numerator += (P[points[i][0]][points[i][1]]*s[i]*((nor[i][0])/(points[i][0] - ind[0])*dx + (nor[i][1])/(points[i][1] - ind[1])*dx));
            denominator += (s[i]*((nor[i][0])/(points[i][0] - ind[0])*dx +(nor[i][1])/(points[i][1] - ind[1])*dx));
        }
        cout<<"n: "<<numerator<< " d: "<<denominator<<endl;
    }   
    P[ind[0]][ind[1]] = numerator/denominator;
}

// void fvm(vector<int> ind){ 
//     double P1; 
//     vector<double> s;
//     vector<vector<int> > points; 
//     vector<vector<double> > center, nor;
//     vector<double> delt;
    
//     vector<int> temp;
//     if(ind[1]==1){
//         temp.push_back(ind[0]);
//         temp.push_back(ind[1]-1);
//     }else{
//         temp.push_back(ind[0]);
//         temp.push_back(ind[1]-2);
//     }  
//     points.push_back(temp);
//     s.push_back((ind[0]+1) - (ind[0]-1));
//     temp.clear();

//     if(ind[0]==1){
//         temp.push_back(ind[0]-1);
//         temp.push_back(ind[1]);
//     }else{
//         temp.push_back(ind[0]-2);
//         temp.push_back(ind[1]);
//     }
//     points.push_back(temp);
//     s.push_back((ind[1]+1) - (ind[1]-1));
//     temp.clear();
//     if(ind[1]==(n-2)){
//         temp.push_back(ind[0]);
//         temp.push_back(ind[1]+1);
//     }else{
//         temp.push_back(ind[0]);
//         temp.push_back(ind[1]+2);
//     }
//     points.push_back(temp);
//     s.push_back((ind[0]+1) - (ind[0]-1));
//     temp.clear();
//     if(ind[0]==(n-2)){
//         temp.push_back(ind[0]+1);
//         temp.push_back(ind[1]);
//     }else{
//         temp.push_back(ind[0]+2);
//         temp.push_back(ind[1]);
//     }
//     s.push_back((ind[1]+1) - (ind[1]-1));
//     points.push_back(temp);
//     temp.clear();

//     for(int i=0; i<points.size(); i++){
//         cout<<points[i][0]<<" "<<points[i][1]<<endl;
//         cout<<P[points[i][0]][points[i][1]]<<endl;
//     }

//     for(int i=0; i<points.size(); i++){
        
//         vector<double> temp;
//         if(i == 0){
//             temp.push_back(dx*(points[0][0] +points[points.size()-1][0])/2);
//             temp.push_back(dx*(points[0][1] +points[points.size()-1][1])/2);
//         }else{
//             temp.push_back(dx*(points[i][0] +points[i-1][0])/2);
//             temp.push_back(dx*(points[i][1] +points[i-1][1])/2);
//         }

//         center.push_back(temp);
//         temp.clear();
//         cout<<"Center"<<endl;
//         cout<<center[i][0]<<" "<<center[i][1]<<endl;
//     }

//     for(int i=0; i<center.size(); i++){
//         vector<double> temp;
//         if(i == center.size()-1){
//             temp = normal(row_return(center, i), row_return(center, 0));
//         }
//         else{
//             temp = normal(row_return(center, i), row_return(center, i+1));
//         }
//         nor.push_back(temp);
//         cout << "***";
//         cout<<nor[i][0]<<" "<<nor[i][1]<<endl;
//         temp.clear();
//         cout << "***";
//     }
//     //length 
//     // for(int i=0; i<center.size(); i++){
//     //     if(i == center.size()-1){
//     //         s.push_back(pow((pow((center[i][0] - center[0][0]), 2) + pow((center[i][1] - center[0][1]), 2)) , 0.5));
//     //     }else{
//     //         s.push_back(pow((pow((center[i][0] - center[i+1][0]), 2) + pow((center[i][1] - center[i+1][1]), 2)) , 0.5));
//     //     }
//     //     cout<<i<< " "<<s[i]<<endl;
//     // } 
     
//     double numerator = 0; 
//     double denominator = 0; 
//     for(int i=0; i<points.size(); i++){
//         cout<<"x: "<<points[i][0]<< " y: "<<points[i][1]<<endl;
//         if(((points[i][0] - ind[0]) == 0) and ((points[i][1] - ind[1]) == 0)){
//             continue;
//         }else if((points[i][0] - ind[0]) == 0){
//             numerator += (P[points[i][0]][points[i][1]]*nor[i][1]*s[i])/((points[i][1] - ind[1])*dx);
//             cout<<"1: "<<P[points[i][0]][points[i][1]]<<" "<<nor[i][1]<<" "<<s[i]<<" "<<(points[i][1] - ind[1])*dx<<endl;
//             denominator += (s[i]*nor[i][1])/((points[i][1] - ind[1])*dx);
//         }else if((points[i][1] - ind[1]) == 0){
//             cout<<"2: "<<P[points[i][0]][points[i][1]]<<" "<<nor[i][0]<<" "<<s[i]<<" "<<(points[i][0] - ind[0])*dx<<endl;
//             numerator += (P[points[i][0]][points[i][1]]*nor[i][0]*s[i])/((points[i][0] - ind[0])*dx);
//             denominator += (s[i]*nor[i][0])/((points[i][0] - ind[0])*dx);
//         }else{
//             cout<<"3: "<<P[points[i][0]][points[i][1]]<<" "<<nor[i][0]<<" "<<s[i]<<" "<<(points[i][0] - ind[0])*dx<<endl;
//             cout<<"4: "<<P[points[i][0]][points[i][1]]<<" "<<nor[i][1]<<" "<<points[i][1]<<" "<<ind[1]*dx<<endl;
//             numerator += (P[points[i][0]][points[i][1]]*s[i]*((nor[i][0])/(points[i][0] - ind[0])*dx + (nor[i][1])/(points[i][1] - ind[1])*dx));
//             denominator += (s[i]*((nor[i][0])/(points[i][0] - ind[0])*dx +(nor[i][1])/(points[i][1] - ind[1])*dx));
//         }
//         cout<<"n: "<<numerator<< " d: "<<denominator<<endl;
//     }   
//     P[ind[0]][ind[1]] = numerator/denominator;
// }

void center_average(vector<int> ind){
    double P1; 
    // vector<double> s;
    vector<vector<int> > points; 
    // vector<vector<double> > center, nor;
    double sum=0;
    
    vector<int> temp;
    temp.push_back(ind[0]-1);
    temp.push_back(ind[1]-1);
    points.push_back(temp);
    temp.clear();
    temp.push_back(ind[0]-1);
    temp.push_back(ind[1]+1);
    points.push_back(temp);
    temp.clear();
    temp.push_back(ind[0]+1);
    temp.push_back(ind[1]+1);
    points.push_back(temp);
    temp.clear();
    temp.push_back(ind[0]+1);
    temp.push_back(ind[1]-1);
    points.push_back(temp);
    temp.clear();

    for(int i=0; i<points.size(); i++){
        sum += P[points[i][0]][points[i][1]];
    }
    P[ind[0]][ind[1]] = sum/points.size();

}

int main(){


    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(j==(n-1)){
                P[i][j]=1;
            }else{
                P[i][j]=0;
            }
            P_copy[i][j]=0;
        }
    }
    // for(int i=0; i<n; i++){
    //     for(int j=0; j<n; j++){
    //         cout<<P[j][i]<<" ";
    //     }
    //     cout<<endl;
    // }


    vector<int> ind;
    // for(int i =1; i<n-1; i++){
    //         for(int j =1; j<n-1; j++){
    //             ind.push_back(i);
    //             ind.push_back(j);
    //             center_average(ind);
    //             ind.clear();
    //         }
    // }
    // for(int i=1; i<n-1; i+=2){
    //         for(int j=1; j<n-1; j+=2){
    //             ind.push_back(i);
    //             ind.push_back(j);
    //             fvm(ind);
    //             ind.clear();
    //             // cout<<i<<" "<<j<<endl;
    //         }
    // }
    // for(int i=0; i<n; i++){
    //     for(int j=0; j<n; j++){
    //         cout<<P[j][i]<<" ";
    //     }
    //     cout<<endl;
    // }
    // ind.push_back(1);
    // ind.push_back(1);
    // fvm(ind);
    // ind.clear();
        // for(int i=0; i<n; i++){
        //     for(int j=0; j<n; j++){ 
        //         P_copy[i][j] = P[i][j];
        //     }
        // }
        for(int i =1; i<n-1; i+=2){
                for(int j =1; j<n-1; j+=2){
                    ind.push_back(i);
                    ind.push_back(j);
                    center_average(ind);
                    ind.clear();
                    // cout<<i<<" "<<j<<endl;
                    
                }
        }
        for(int i =2; i<n-1; i+=2){
            for(int j =2; j<n-1; j+=2){
                ind.push_back(i);
                ind.push_back(j);
                fvm(ind);
                ind.clear();
                // cout<<i<<" "<<j<<endl;
                
            }
        }
    while((abs(find_max(P_copy)-find_max(P))>e) ){
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){ 
                P_copy[i][j] = P[i][j];
            }
        }
        for(int i =1; i<n-1; i+=2){
                for(int j =1; j<n-1; j+=2){
                    ind.push_back(i);
                    ind.push_back(j);
                    center_average(ind);
                    ind.clear();
                    // cout<<i<<" "<<j<<endl;
                    
                }
        }
        for(int i =2; i<n-1; i+=2){
            for(int j =2; j<n-1; j+=2){
                ind.push_back(i);
                ind.push_back(j);
                fvm(ind);
                ind.clear();
                // cout<<i<<" "<<j<<endl;
                
            }
        }
    }

    for(int i=1; i<(n-1); i++){
        for(int j=1; j<(n-1); j++){
            if(P[i][j]==0){
                ind.push_back(i);
                ind.push_back(j);
                center_average(ind);
                ind.clear();
            }
        }
        // cout<<endl;
    }

    // while((abs(find_max(P_copy)-find_max(P))>e) ){
    //     for(int i=0; i<n; i++){
    //         for(int j=0; j<n; j++){ 
    //             P_copy[i][j] = P[i][j];
    //         }
    //     }
    //     for(int i =1; i<n-1; i++){
    //         for(int j =1; j<n-1; j++){
    //             ind.push_back(i);
    //             ind.push_back(j);
    //             fvm(ind);
    //             ind.clear();
    //         }
    //     }
    //     cout<<find_max(P_copy)<<" "<<find_max(P)<<endl;
    // }

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(i % 2 == 1 or j % 2 == 1){
                continue;
            }else{
            cout<<P[j][i]<<" ";
            }
        }
        cout<<endl;

    }

    return 0;
}