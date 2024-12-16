#include <bits/stdc++.h>
using namespace std;
void check(int array[]){
    cout<<array[1]<<endl;
}
int main(){
    int cnt = 0;
    int array[12];
    for(int i=0; i<12; i++){
        cnt += 30;
        array[i] = cnt;
        cout<<"degrees: "<<cnt<<" "<<"cos: "<<cos(cnt*(M_PI/180))<<" sin: "<<sin(cnt*(M_PI/180))<<endl;
    }
    check(array);
    return 0;
}