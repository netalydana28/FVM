#include <bits/stdc++.h>  
#include <string>
using namespace std; 
vector<double> split_double(string line, char delimeter){
    string my_str = "";
    vector<double> v;
    for(int i=0; i<line.size(); i++){
        if(((line[i] == delimeter) or (i==line.size()-1))){
            if(line[i+1] == delimeter){
                continue;
            }
            if(i==line.size()-1){
                my_str+=line[i];  
            }
            char* temp = new char[my_str.length()+1];
            strcpy(temp, my_str.c_str());
            double to_insert = atof(temp);
            v.push_back(to_insert);
            my_str = "";
        }else{
            my_str+=line[i];  
        }
    }
    return v;
}

void check_error(string right, string check){
    string right_Text, check_Text, first = "out" + right + ".dat", second = "out" + check + ".dat";
    ifstream rightFile(first.c_str());
    ifstream checkFile(second.c_str());
    vector<double> v_right, v_check; 
    int cnt = 1;
    double max_diff_P = 0.0, max_diff_C = 0.0, max_diff_U = 0.0, max_diff_V = 0.0, 
        node_num_P = 0.0, node_num_U = 0.0, node_num_V = 0.0, node_num_C = 0.0; 
    while (getline(rightFile, right_Text)) {
        v_right = split_double(right_Text, ' ');
        v_right.erase(v_right.begin());
        v_right.erase(v_right.begin());
        getline(checkFile, check_Text);
        v_check = split_double(check_Text, ' ');
        v_check.erase(v_check.begin());
        v_check.erase(v_check.begin());
        if(fabs(v_right[0] - v_check[0]) > max_diff_P){
            node_num_P = cnt; 
            max_diff_P = fabs(v_right[0] - v_check[0]);
        }
        if(fabs(v_right[1] - v_check[1]) > max_diff_U){
            node_num_U = cnt; 
            max_diff_U = fabs(v_right[1] - v_check[1]);
        }
        if(fabs(v_right[2] - v_check[2]) > max_diff_V){
            node_num_V = cnt; 
            max_diff_V = fabs(v_right[2] - v_check[2]);
        }
        if(fabs(v_right[3] - v_check[3]) > max_diff_C){
            node_num_C = cnt; 
            max_diff_C = fabs(v_right[3] - v_check[3]);
        }
        v_right.clear();
        v_check.clear();
        cnt++;
    }
    rightFile.close();
    checkFile.close();
    printf("%s iterations:\n node num %.0f, max diff P: %.6f\n  node num %.0f, max diff U: %.6f\n   node num %.0f, max diff V: %.6f\n    node num %.0f, max diff C: %.6f\n", right.c_str(), node_num_P, max_diff_P, node_num_U, max_diff_U, node_num_V, max_diff_V, node_num_C, max_diff_C);
}
int main(){
    check_error("500K", "500000");
    check_error("1M", "1000000");
    check_error("1M500K", "1500000");
    check_error("2M", "2000000");
    check_error("2M500K", "2500000");

    return 0;
}