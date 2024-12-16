    #include <bits/stdc++.h>  
using namespace std;   
  
vector<vector<double> > coor;

vector<double> row_return(vector<vector<double> > v, int n){
    vector<double> row; 
    for(int i; i<v[n].size(); i++){
        row.push_back(v[n][i]);
    }
    return row;
}

// vector<double> normal(vector<double> l1, vector<double> l2, string direction){ 
//     vector <double> ans; 
//     cout<<"   "<<l1[0]<<" "<<l1[0]<< " "<<l2[0]<<" "<<l2[0]<<endl;
//     if (direction == "+"){ 
//         ans.push_back(-((l1[1] - l2[1]) + (l1[1] - l2[1]))); 
//         ans.push_back(-((l2[0] - l1[0]) + (l2[0] - l1[0])));  
//     }else if (direction == "-"){ 
//         ans.push_back((l1[1] - l2[1]) + (l1[1] - l2[1])); 
//         ans.push_back((l2[0] - l1[0]) + (l2[0] - l1[0]));  
//     } 
//     double len = pow(pow(ans[0], 2) + pow(ans[1], 2), 0.5); 
//     for(int i=0; i<ans.size(); i++){
//         ans[i] /= len;
//     }
//     //cout<<ans[0]<<" "<<ans[1]<<endl; 
//     return ans; 
// } 

vector<double> normal(vector<double> l1, vector<double> l2){ 
    vector <double> ans; 
   
    //cout << l1[0] << " " << l1[1] << " " << l2[0] << " " << l2[1] << endl;
    ans.push_back(l1[1] - l2[1]); 
    ans.push_back(l2[0] - l1[0]);  

    double len = pow(pow(ans[0], 2) + pow(ans[1], 2), 0.5); 
    for(int i=0; i<ans.size(); i++){
        ans[i] /= len;
    }
    return ans; 
}
 

// vector<double> get_vertex(vector<vector<double> > points){
//     //finding equation 
//     double m1, m2, b1, b2; 
//     m1 = (points[0][1] - points[2][1])/(points[0][0] - points[2][0]);
//     m2 = (points[1][1] - points[1][1])/(points[1][0] - points[3][0]); 
//     b1 = points[0][1] - m1 * points[0][0]; 
//     b2 = points[1][1] - m1 * points[1][0]; 
//     vector<double> intersection; 
//     intersection.push_back((b2-b1)/(m1-m2));
//     intersection.push_back(m1*intersection[0] - b1);
//     return intersection;
// }

 
double fvm(vector<vector<double> > points, vector<double> P){ 
    double P1; 
    vector<double>  s;
    vector<vector<double> > center, n;
    //cout<<sizeof(ind)/sizeof(ind[0])<<"N";
//     for(int i=0; i<ind.size(); i++){
//         if(i==ind.size()-1){
//             vertex.push_back(coor[ind[i]][0]);
//             vertex.push_back(coor[ind[i]][1]);
//         }else{
//             points.push_back(coor[ind[i]]);
// //            cout<<"Points"<<endl;
// //            cout<<points[i][0]<<" "<<points[i][1]<<endl;
//         }
//     }
    // x1 = points[ind2].first - points[ind4].first;  
    // y1 = points[ind3].second - points[ind5].second; 
    //center determination
    //cout<<points.size()<<"!";
    // for(int i=0; i<points.size(); i++){
    //     //cout<<i;
    //     vector<double> temp;
    //     if(i == points.size()-1){
    //         temp.push_back((points[0][0] +points[i][0])/2);
    //         temp.push_back((points[0][1] +points[i][1])/2);
    //     }else{
    //         temp.push_back((points[i][0] +points[i+1][0])/2);
    //         temp.push_back((points[i][1] +points[i+1][1])/2);
    //     }
    for(int i=0; i<points.size(); i++){
        //cout<<i;
        vector<double> temp;
        if(i == 0){
            temp.push_back((points[0][0] +points[points.size()-1][0])/2);
            temp.push_back((points[0][1] +points[points.size()-1][1])/2);
        }else{
            temp.push_back((points[i][0] +points[i-1][0])/2);
            temp.push_back((points[i][1] +points[i-1][1])/2);
        }

        //cout<<points[i][0]<<" "<<points[i][1]<<endl;
        center.push_back(temp);
        temp.clear();
//        cout<<"Center"<<endl;
//        cout<<center[i][0]<<" "<<center[i][1]<<endl;
    }
    //normal vector determination
    // for(int i=0; i<center.size(); i++){
    //     vector<double> temp;
//     if(i == center.size()-1){
    //         // cout<<"here1"<<center[i][0]<<" "<<center[0][0]<<endl;
    //         // cout<<"here2"<<center[i][1]<<" "<<center[0][1]<<endl;
    //         if((center[i][0]<vertex[0] and center[0][0]<vertex[0]) or (center[i][1]<vertex[1] and center[0][1]<vertex[1])){
    //             temp = normal(row_return(center, i), row_return(center, 0), "-");
    //         }else{
    //             temp = normal(row_return(center, i), row_return(center, 0), "+");
    //         }
    //     }else{
    //         // cout<<"here1"<<center[i][0]<<" "<<center[i+1][0]<<endl;
    //         // cout<<"here2"<<center[i][1]<<" "<<center[i+1][1]<<endl;
    //         if((center[i][0]<vertex[0] and center[i+1][0]<vertex[0]) or (center[i][1]<vertex[1] and center[i+1][1]<vertex[1])){
    //             temp = normal(row_return(center, i), row_return(center, i+1), "-");
    //         }else{
    //             temp = normal(row_return(center, i), row_return(center, i+1), "+");
    //         }
    //     }
    //     n.push_back(temp);
    //     cout<<n[i][0]<<" "<<n[i][1]<<endl;
    //     temp.clear();
    // }
    for(int i=0; i<center.size(); i++){
        vector<double> temp;
        if(i == center.size()-1){
            temp = normal(row_return(center, i), row_return(center, 0));
        }
        else{
            temp = normal(row_return(center, i), row_return(center, i+1));
        }
        n.push_back(temp);
//        cout << "***";
//        cout<<n[i][0]<<" "<<n[i][1]<<endl;
        temp.clear();
//        cout << "***";
    }
    //length 
    for(int i=0; i<center.size(); i++){
        if(i == center.size()-1){
            s.push_back(pow((pow((center[i][0] - center[0][0]), 2) + pow((center[i][1] - center[0][1]), 2)) , 0.5));
        }else{
            s.push_back(pow((pow((center[i][0] - center[i+1][0]), 2) + pow((center[i][1] - center[i+1][1]), 2)) , 0.5));
        }
//        cout<<i<< " "<<s[i]<<endl;
    } 
     
    double numerator = 0; 
    double denominator = 0; 
    for(int i=0; i<points.size(); i++){
        if(((points[i][0] - points[4][0]) == 0) and ((points[i][1] - points[4][1]) == 0)){
            continue;
        }else if((points[i][0] - points[4][0]) == 0){
            numerator += (P[i]*n[i][1]*s[i])/(points[i][1] - points[4][1]);
//            cout<<"1: "<<P[i]<<" "<<n[i][1]<<" "<<points[i][1]<<" "<<vertex[1]<<endl;
            denominator += (s[i]*n[i][1])/(points[i][1] - points[4][1]);
        }else if((points[i][1] - points[4][1]) == 0){
//            cout<<"1: "<<P[i]<<" "<<n[i][0]<<" "<<points[i][0]<<" "<<vertex[0]<<endl;
            numerator += (P[i]*n[i][0]*s[i])/(points[i][0] - points[4][0]);
            denominator += (s[i]*n[i][0])/(points[i][0] - points[4][0]);
        }else{
//            cout<<"1: "<<P[i]<<" "<<n[i][0]<<" "<<points[i][0]<<" "<<vertex[0]<<endl;
//            cout<<"2: "<<P[i]<<" "<<n[i][1]<<" "<<points[i][1]<<" "<<vertex[1]<<endl;
            numerator += (P[i]*s[i]*((n[i][0])/(points[i][0] - points[4][0]) + (n[i][1])/(points[i][1] - points[4][1])));
            denominator += (s[i]*((n[i][0])/(points[i][0] - points[4][0]) +(n[i][1])/(points[i][1] - points[4][1])));
        }
//        cout<<"n: "<<numerator<< " d: "<<denominator<<endl;
    }   
    return numerator/denominator;
}  
   
int main()  
{   
    // Initializing the vector of vectors
  
    int n;  
    cin>>n;  
    double x, y;  
    x = 0;  
    y = 0;  
    double dx = 0.1;
    double dy = 0.1; 
  
    // Inserting elements into vector
    for (double i = 0; i < 1; i+=0.1) {
        vector<double> v1;
        for (double j = 0; j < 1; j+=0.1) {
            v1.push_back(j);
            v1.push_back(i);
            coor.push_back(v1);
            v1.clear();
        }
        
    }
    cout << coor.size() << "******************";
    for (int i = 0; i < coor.size(); i++) {
        for (int j = 0; j < coor[i].size(); j++)
            cout << coor[i][j] << " ";
        cout << endl;
    }
    cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&";

    vector<int> ind; //= {1, 10, 21, 12};
    ind.push_back(1);
    ind.push_back(10);
    ind.push_back(21);
    ind.push_back(12);
    ind.push_back(11);
    vector<double> P;
    P.push_back(15);
    P.push_back(10);
    P.push_back(20);
    P.push_back(25);
    
    //double P1; 
    //P1 = fvm(ind, P); 
//    cout<<P1<<endl;
    
    int n0 = 1 / dy;
    cout << n0<<"(((";
    int m0 = 1 / dx;
    double matrix[n0][m0];
    for (int i = 0; i < n0; i++) {
        for (int j = 0; j < m0; j++){
            if(j == 0){
                matrix[i][j] = 1;
            }else{
                matrix[i][j] = 0;
            }
        }
    }
      
    for (int i = 1; i < n0-1; i++) {
        for (int j = 1; j < m0-1; j++){
            vector<vector<double> > points0;
            vector<double> point1 = {i*0.1, (j+1)*0.1};
            cout << i*0.1 << " " << (j-1)*0.1;
            vector<double> point2 = {(i-1)*0.1, j*0.1};
            vector<double> point3 = {i*0.1, (j-1)*0.1};
            vector<double> point4 = {(i+1)*0.1, j*0.1};
            vector<double> point5 = {i*0.1, j*0.1};
            points0 = {point3, point2, point1, point4, point5};
            vector<double> pressures = {matrix[i][j-1], matrix[i-1][j], matrix[i][j+1], matrix[i+1][j]};
            matrix[i][j] = fvm(points0, pressures);
            cout << "=========================" << matrix[i][j] << " " << i << " " << j << endl;
        }
    }
    cout << n0<< " " << m0<< endl;
    for (int i = 0; i < n0; i++) {
        for (int j = 0; j < m0; j++){
            cout << matrix[i][j] << " "; 
        }
        cout << endl;
    }
      
    return 0;  
}    #include <bits/stdc++.h>  
using namespace std;   
  
vector<vector<double> > coor;

vector<double> row_return(vector<vector<double> > v, int n){
    vector<double> row; 
    for(int i; i<v[n].size(); i++){
        row.push_back(v[n][i]);
    }
    return row;
}

// vector<double> normal(vector<double> l1, vector<double> l2, string direction){ 
//     vector <double> ans; 
//     cout<<"   "<<l1[0]<<" "<<l1[0]<< " "<<l2[0]<<" "<<l2[0]<<endl;
//     if (direction == "+"){ 
//         ans.push_back(-((l1[1] - l2[1]) + (l1[1] - l2[1]))); 
//         ans.push_back(-((l2[0] - l1[0]) + (l2[0] - l1[0])));  
//     }else if (direction == "-"){ 
//         ans.push_back((l1[1] - l2[1]) + (l1[1] - l2[1])); 
//         ans.push_back((l2[0] - l1[0]) + (l2[0] - l1[0]));  
//     } 
//     double len = pow(pow(ans[0], 2) + pow(ans[1], 2), 0.5); 
//     for(int i=0; i<ans.size(); i++){
//         ans[i] /= len;
//     }
//     //cout<<ans[0]<<" "<<ans[1]<<endl; 
//     return ans; 
// } 

vector<double> normal(vector<double> l1, vector<double> l2){ 
    vector <double> ans; 
   
    //cout << l1[0] << " " << l1[1] << " " << l2[0] << " " << l2[1] << endl;
    ans.push_back(l1[1] - l2[1]); 
    ans.push_back(l2[0] - l1[0]);  

    double len = pow(pow(ans[0], 2) + pow(ans[1], 2), 0.5); 
    for(int i=0; i<ans.size(); i++){
        ans[i] /= len;
    }
    return ans; 
}
 

// vector<double> get_vertex(vector<vector<double> > points){
//     //finding equation 
//     double m1, m2, b1, b2; 
//     m1 = (points[0][1] - points[2][1])/(points[0][0] - points[2][0]);
//     m2 = (points[1][1] - points[1][1])/(points[1][0] - points[3][0]); 
//     b1 = points[0][1] - m1 * points[0][0]; 
//     b2 = points[1][1] - m1 * points[1][0]; 
//     vector<double> intersection; 
//     intersection.push_back((b2-b1)/(m1-m2));
//     intersection.push_back(m1*intersection[0] - b1);
//     return intersection;
// }

 
double fvm(vector<vector<double> > points, vector<double> P){ 
    double P1; 
    vector<double>  s;
    vector<vector<double> > center, n;
    //cout<<sizeof(ind)/sizeof(ind[0])<<"N";
//     for(int i=0; i<ind.size(); i++){
//         if(i==ind.size()-1){
//             vertex.push_back(coor[ind[i]][0]);
//             vertex.push_back(coor[ind[i]][1]);
//         }else{
//             points.push_back(coor[ind[i]]);
// //            cout<<"Points"<<endl;
// //            cout<<points[i][0]<<" "<<points[i][1]<<endl;
//         }
//     }
    // x1 = points[ind2].first - points[ind4].first;  
    // y1 = points[ind3].second - points[ind5].second; 
    //center determination
    //cout<<points.size()<<"!";
    // for(int i=0; i<points.size(); i++){
    //     //cout<<i;
    //     vector<double> temp;
    //     if(i == points.size()-1){
    //         temp.push_back((points[0][0] +points[i][0])/2);
    //         temp.push_back((points[0][1] +points[i][1])/2);
    //     }else{
    //         temp.push_back((points[i][0] +points[i+1][0])/2);
    //         temp.push_back((points[i][1] +points[i+1][1])/2);
    //     }
    for(int i=0; i<points.size(); i++){
        //cout<<i;
        vector<double> temp;
        if(i == 0){
            temp.push_back((points[0][0] +points[points.size()-1][0])/2);
            temp.push_back((points[0][1] +points[points.size()-1][1])/2);
        }else{
            temp.push_back((points[i][0] +points[i-1][0])/2);
            temp.push_back((points[i][1] +points[i-1][1])/2);
        }

        //cout<<points[i][0]<<" "<<points[i][1]<<endl;
        center.push_back(temp);
        temp.clear();
//        cout<<"Center"<<endl;
//        cout<<center[i][0]<<" "<<center[i][1]<<endl;
    }
    //normal vector determination
    // for(int i=0; i<center.size(); i++){
    //     vector<double> temp;
//     if(i == center.size()-1){
    //         // cout<<"here1"<<center[i][0]<<" "<<center[0][0]<<endl;
    //         // cout<<"here2"<<center[i][1]<<" "<<center[0][1]<<endl;
    //         if((center[i][0]<vertex[0] and center[0][0]<vertex[0]) or (center[i][1]<vertex[1] and center[0][1]<vertex[1])){
    //             temp = normal(row_return(center, i), row_return(center, 0), "-");
    //         }else{
    //             temp = normal(row_return(center, i), row_return(center, 0), "+");
    //         }
    //     }else{
    //         // cout<<"here1"<<center[i][0]<<" "<<center[i+1][0]<<endl;
    //         // cout<<"here2"<<center[i][1]<<" "<<center[i+1][1]<<endl;
    //         if((center[i][0]<vertex[0] and center[i+1][0]<vertex[0]) or (center[i][1]<vertex[1] and center[i+1][1]<vertex[1])){
    //             temp = normal(row_return(center, i), row_return(center, i+1), "-");
    //         }else{
    //             temp = normal(row_return(center, i), row_return(center, i+1), "+");
    //         }
    //     }
    //     n.push_back(temp);
    //     cout<<n[i][0]<<" "<<n[i][1]<<endl;
    //     temp.clear();
    // }
    for(int i=0; i<center.size(); i++){
        vector<double> temp;
        if(i == center.size()-1){
            temp = normal(row_return(center, i), row_return(center, 0));
        }
        else{
            temp = normal(row_return(center, i), row_return(center, i+1));
        }
        n.push_back(temp);
//        cout << "***";
//        cout<<n[i][0]<<" "<<n[i][1]<<endl;
        temp.clear();
//        cout << "***";
    }
    //length 
    for(int i=0; i<center.size(); i++){
        if(i == center.size()-1){
            s.push_back(pow((pow((center[i][0] - center[0][0]), 2) + pow((center[i][1] - center[0][1]), 2)) , 0.5));
        }else{
            s.push_back(pow((pow((center[i][0] - center[i+1][0]), 2) + pow((center[i][1] - center[i+1][1]), 2)) , 0.5));
        }
//        cout<<i<< " "<<s[i]<<endl;
    } 
     
    double numerator = 0; 
    double denominator = 0; 
    for(int i=0; i<points.size(); i++){
        if(((points[i][0] - points[4][0]) == 0) and ((points[i][1] - points[4][1]) == 0)){
            continue;
        }else if((points[i][0] - points[4][0]) == 0){
            numerator += (P[i]*n[i][1]*s[i])/(points[i][1] - points[4][1]);
//            cout<<"1: "<<P[i]<<" "<<n[i][1]<<" "<<points[i][1]<<" "<<vertex[1]<<endl;
            denominator += (s[i]*n[i][1])/(points[i][1] - points[4][1]);
        }else if((points[i][1] - points[4][1]) == 0){
//            cout<<"1: "<<P[i]<<" "<<n[i][0]<<" "<<points[i][0]<<" "<<vertex[0]<<endl;
            numerator += (P[i]*n[i][0]*s[i])/(points[i][0] - points[4][0]);
            denominator += (s[i]*n[i][0])/(points[i][0] - points[4][0]);
        }else{
//            cout<<"1: "<<P[i]<<" "<<n[i][0]<<" "<<points[i][0]<<" "<<vertex[0]<<endl;
//            cout<<"2: "<<P[i]<<" "<<n[i][1]<<" "<<points[i][1]<<" "<<vertex[1]<<endl;
            numerator += (P[i]*s[i]*((n[i][0])/(points[i][0] - points[4][0]) + (n[i][1])/(points[i][1] - points[4][1])));
            denominator += (s[i]*((n[i][0])/(points[i][0] - points[4][0]) +(n[i][1])/(points[i][1] - points[4][1])));
        }
//        cout<<"n: "<<numerator<< " d: "<<denominator<<endl;
    }   
    return numerator/denominator;
}  
   
int main()  
{   
    // Initializing the vector of vectors
  
    int n;  
    cin>>n;  
    double x, y;  
    x = 0;  
    y = 0;  
    double dx = 0.1;
    double dy = 0.1; 
  
    // Inserting elements into vector
    for (double i = 0; i < 1; i+=0.1) {
        vector<double> v1;
        for (double j = 0; j < 1; j+=0.1) {
            v1.push_back(j);
            v1.push_back(i);
            coor.push_back(v1);
            v1.clear();
        }
        
    }
    cout << coor.size() << "******************";
    for (int i = 0; i < coor.size(); i++) {
        for (int j = 0; j < coor[i].size(); j++)
            cout << coor[i][j] << " ";
        cout << endl;
    }
    cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&";

    vector<int> ind; //= {1, 10, 21, 12};
    ind.push_back(1);
    ind.push_back(10);
    ind.push_back(21);
    ind.push_back(12);
    ind.push_back(11);
    vector<double> P;
    P.push_back(15);
    P.push_back(10);
    P.push_back(20);
    P.push_back(25);
    
    //double P1; 
    //P1 = fvm(ind, P); 
//    cout<<P1<<endl;
    
    int n0 = 1 / dy;
    cout << n0<<"(((";
    int m0 = 1 / dx;
    double matrix[n0][m0];
    for (int i = 0; i < n0; i++) {
        for (int j = 0; j < m0; j++){
            if(j == 0){
                matrix[i][j] = 1;
            }else{
                matrix[i][j] = 0;
            }
        }
    }
      
    for (int i = 1; i < n0-1; i++) {
        for (int j = 1; j < m0-1; j++){
            vector<vector<double> > points0;
            vector<double> point1 = {i*0.1, (j+1)*0.1};
            cout << i*0.1 << " " << (j-1)*0.1;
            vector<double> point2 = {(i-1)*0.1, j*0.1};
            vector<double> point3 = {i*0.1, (j-1)*0.1};
            vector<double> point4 = {(i+1)*0.1, j*0.1};
            vector<double> point5 = {i*0.1, j*0.1};
            points0 = {point3, point2, point1, point4, point5};
            vector<double> pressures = {matrix[i][j-1], matrix[i-1][j], matrix[i][j+1], matrix[i+1][j]};
            matrix[i][j] = fvm(points0, pressures);
            cout << "=========================" << matrix[i][j] << " " << i << " " << j << endl;
        }
    }
    cout << n0<< " " << m0<< endl;
    for (int i = 0; i < n0; i++) {
        for (int j = 0; j < m0; j++){
            cout << matrix[i][j] << " "; 
        }
        cout << endl;
    }
      
    return 0;  
}