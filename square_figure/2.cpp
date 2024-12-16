#include <bits/stdc++.h>  
using namespace std;   
  
map <int, pair<double, points; 

pair<double, double> normal(pair<double, double> l1, pair<double, double> l2, string direction){ 
    double ans[2];; 
    if (direction == "+"){ 
        ans[0] = l1.second - l2.second; 
        ans[1] = l2.first - l1.first;  
    }else if (direction == "-"){ 
        ans[0] = -(l1.second - l2.second); 
        ans[1] = -(l2.first - l1.first);  
    } 
    //cout<<ans[0]<<" "<<ans[1]<<endl; 
    return make_pair(ans[0], ans[1]); 
} 
 
// double* fvm(double l1[2], double l2[2], double l3[2], double l4[2],, double l5[2], P2, P3, P4, P5){ 
     
// } 
 
double length(pair<double, double> l3, pair<double, double> l2, pair<double, double> l5){ 
    pair<double, double> center23; 
    center23.first = (l3.first + l2.first) / 2; 
    center23.second = (l3.second + l2.second) / 2; 
    pair<double, double> center25; 
    center25.first = (l5.first + l2.first) / 2; 
    center25.second = (l5.second + l2.second) / 2; 
    double len = pow((pow((center23.first - center25.first), 2) + pow((center23.second - center25.second), 2)) , 0.5); 
    return len; 
} 
 
void fvm(int ind1, int ind2, double P2, int ind3, double P3, int ind4, double P4, int ind5, double P5){ 
        double P1; 
        pair<double, double> n2, n3, n4, n5; 
        // x1 = points[ind2].first - points[ind4].first;  
        // y1 = points[ind3].second - points[ind5].second;  
        n2 = normal(points[ind2], points[ind4], "-"); 
        n3 = normal(points[ind3], points[ind5], "+"); 
        n4 = normal(points[ind3], points[ind5], "-"); 
        n5 = normal(points[ind2], points[ind4], "+"); 
        //return P1; 
        double s21 = length(points[ind3], points[ind2], points[ind5]); 
        double s31 = length(points[ind4], points[ind3], points[ind2]); 
        double s41 = length(points[ind5], points[ind4], points[ind3]); 
        double s51 = length(points[ind2], points[ind5], points[ind4]); 
        cout << s21 << ", " << s31 << ", " << s41 << ", " << s51<<endl; 
        cout << n2.first << ", " << n2.second<<endl; 
        cout << n3.first << ", " << n3.second<<endl; 
        cout << n4.first << ", " << n4.second<<endl; 
        cout << n5.first << ", " << n5.second<<endl; 
        cout << endl; 
        cout << points[ind1].first << ", " << points[ind1].second<<endl; 
        cout << points[ind2].first << ", " << points[ind2].second<<endl; 
        cout << points[ind3].first << ", " << points[ind3].second<<endl; 
        cout << points[ind4].first << ", " << points[ind4].second<<endl; 
        cout << points[ind5].first << ", " << points[ind5].second<<endl; 
         
        P1 = P2*(n2.first/(points[ind2].first - points[ind1].first) + n2.second/(points[ind2].second - points[ind1].second))*s21 + P3*(n3.first/(points[ind3].first - points[ind1].first) + n3.second/(points[ind3].second - points[ind1].second))*s31 + P4*(n4.first/(points[ind4].first - points[ind1].first) + n4.second/(points[ind4].second - points[ind1].second))*s41 + P5*(n5.first/(points[ind5].first - points[ind1].first) + n5.second/(points[ind5].second - points[ind1].second))*s51; 
        cout << P1; 
        P1 /= (s21*n2.first/(points[ind2].first - points[ind1].first) + s21*n2.second/(points[ind2].second - points[ind1].second) + s31*n3.first/(points[ind3].first - points[ind1].first) + s31*n3.second/(points[ind3].second - points[ind1].second) + s41*n4.first/(points[ind4].first - points[ind1].first) + s41*n4.second/(points[ind4].second - points[ind1].second) + s51*n5.first/(points[ind5].first - points[ind1].first) + s51*n5.second/(points[ind5].second - points[ind1].second)); 
        cout << P1; 
} 
 
   
int main()  
{   
    //cout << "978";  
    int n;  
    cin>>n;  
    //map <int, pair<double, double>> points;   
    double x, y;  
    x = 0;  
    y = 0;  
    int cnt =0;  
    for(double i=0; i<n; i+=0.5){  
         for(double j=0; j<n; j+=0.5){  
            //double coor[2] = {x, y};  
            cnt+=1;  
            points[cnt] = make_pair(i, j);   
            //cout<<i<<" "<<j<<endl; 
              
        } 
    } 
      
        //cout<<points[i].first<<", "<<points[i].second<<endl;  
    //}  
 
    //double p1[2] = {-0.5, -0.5};  
    //double p2[2] = {-0.5, 0.5};  
    //pair<double, double> nor = normal(p1, p2, "+"); 
     
     
    //cout<<nor.first<<" "<<nor.second; 
     
    fvm(12, 2, 123, 11, 1243, 22, 876, 13, 876); 
    cout<<points[1].first<<" "<<points[1].second<<endl; 
     
     
    return 0;  
}