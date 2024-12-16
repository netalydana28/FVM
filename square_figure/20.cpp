// Triangle mesh generator optimization, map to array
#include <bits/stdc++.h>  
using namespace std;
const int n = 80;
const int n_side = n+1;
const int node_num=n_side*n_side;
const double step = 5/double(n);
map<int, vector<double> > nodes;
const int center_num=(n_side-1)*(n_side-1)*2;
map<int, vector<double> > centers;
const int face_num = ((n_side-1)*(n_side)*2) + (center_num/2);
map<vector<int> , int> l_cell_to_face;
map<vector<int> , int> l_face_to_cell;
map<vector<int> , int> l_face_to_node;
map<vector<int> , int> node_to_cell;
map<vector<int> , double> sn;
double P_nodes[node_num];
double P_faces[face_num];
double P_centers[center_num];
double P_centers_copy[center_num];
const double e=0.000001;
double gauss_seidel_coeff[center_num][center_num];
double gauss_seidel_b[center_num];
double find_max(double  v[center_num]){
    double max = v[1];
    for(int i=1; i<center_num; i++){
        if(v[i]>max){
            max = v[i];
        }
    }
    return max;
}
vector<double> normal(int face_ind){
    vector <double> ans; 
    vector <int> key; 
    key.push_back(face_ind); 
    key.push_back(1);
    vector<double> n1 = nodes[l_face_to_node[key]];
    key.erase(key.end()-1);
    key.push_back(2);
    vector<double> n2 = nodes[l_face_to_node[key]];
    ans.push_back(n1[1] - n2[1]); 
    ans.push_back(n2[0] - n1[0]);  
    double len = pow(pow(ans[0], 2) + pow(ans[1], 2), 0.5); 
    for(int i=0; i<ans.size(); i++){
        ans[i] /= len;
    }
    return ans; 
}
double area(int face_ind){
    vector <int> key; 
    key.push_back(face_ind); 
    key.push_back(1);
    vector<double> n1 = nodes[l_face_to_node[key]];
    key.erase(key.end()-1);
    key.push_back(2);
    vector<double> n2 = nodes[l_face_to_node[key]];
    return pow((pow((n1[0] - n2[0]), 2) + pow((n1[1] - n2[1]), 2)) , 0.5);
}
bool is_boundary_face(int face_ind){
    vector<int> key1;
    key1.push_back(face_ind);
    key1.push_back(1);
    vector<int> key2; 
    key2.push_back(face_ind);
    key2.push_back(2);
    if(l_face_to_cell[key1] == l_face_to_cell[key2]){
        return true;
    }
    return false;
}
vector<double> face_center(int face_ind){
    vector <int> key; 
    key.push_back(face_ind); 
    key.push_back(1);
    vector<double> n1 = nodes[l_face_to_node[key]];
    key.erase(key.end()-1);
    key.push_back(2);
    vector<double> n2 = nodes[l_face_to_node[key]];
    vector<double> ans; 
    ans.push_back((n1[0] + n2[0])/2);
    ans.push_back((n1[1] + n2[1])/2);
    return ans;
}
double cell_to_vertex_interpolation(int node_ind){
    vector<double> node_coor = nodes[node_ind]; 
    vector<int> key; 
    key.push_back(node_ind); 
    key.push_back(1); 
    register int cell1 = node_to_cell[key]; 
    key.erase(key.end()-1);  
    key.push_back(2); 
    register int cell2 = node_to_cell[key]; 
    key.erase(key.end()-1);  
    key.push_back(3); 
    register int cell3 = node_to_cell[key]; 
    key.erase(key.end()-1);  
    key.push_back(4); 
    register int cell4 = node_to_cell[key]; 
    key.erase(key.end()-1);  
    key.push_back(5); 
    register int cell5 = node_to_cell[key];
    key.erase(key.end()-1);  
    key.push_back(6); 
    register int cell6 = node_to_cell[key];  
    vector<double> center1 = centers[cell1]; 
    vector<double> center2 = centers[cell2]; 
    vector<double> center3 = centers[cell3]; 
    vector<double> center4 = centers[cell4]; 
    vector<double> center5 = centers[cell5]; 
    vector<double> center6 = centers[cell6];
    double P1 = P_centers[cell1-1]; 
    double P2 = P_centers[cell2-1]; 
    double P3 = P_centers[cell3-1]; 
    double P4 = P_centers[cell4-1]; 
    double P5 = P_centers[cell5-1]; 
    double P6 = P_centers[cell6-1];
    double d1 = pow(pow((center1[0] - node_coor[0]), 2) + pow((center1[1] - node_coor[1]), 2), 0.5); 
    double d2 = pow(pow((center2[0] - node_coor[0]), 2) + pow((center2[1] - node_coor[1]), 2), 0.5); 
    double d3 = pow(pow((center3[0] - node_coor[0]), 2) + pow((center3[1] - node_coor[1]), 2), 0.5); 
    double d4 = pow(pow((center4[0] - node_coor[0]), 2) + pow((center4[1] - node_coor[1]), 2), 0.5); 
    double d5 = pow(pow((center5[0] - node_coor[0]), 2) + pow((center5[1] - node_coor[1]), 2), 0.5); 
    double d6 = pow(pow((center6[0] - node_coor[0]), 2) + pow((center6[1] - node_coor[1]), 2), 0.5); 
    return (P1/d1 + P2/d2 + P3/d3 + P4/d4 + P5/d5 + P6/d6)/(1/d1 + 1/d2 + 1/d3 + 1/d4 + 1/d5 + 1/d6); 
}
void gauss_seidel(){
    double to_check = 1;
    double w = 2/(1+sin(M_PI*(1/(double(center_num+1)))));
    cout<<"w: "<<w<<endl;
    while(abs(to_check - find_max(P_centers))>e){
        to_check = find_max(P_centers);
        for(int i=0; i<center_num; i++){
            P_centers_copy[i] = gauss_seidel_b[i];
            for(int j=0; j<center_num; j++){
                if (j == i){
                    continue;
                }
                P_centers_copy[i] = P_centers_copy[i] - (gauss_seidel_coeff[i][j] * P_centers[j]);
                P_centers[i] = (1* (P_centers_copy[i]/gauss_seidel_coeff[i][i]))+ ((1-1)*P_centers[i]);
            }
        }
    }
}
void fvm(){
    for(int i=1; i<center_num+1; i++){
        vector<double> gDiff;
        vector<int> key;
        key.push_back(i);
        double ac=0;
        double bc=0;
        vector<double> center_coor = centers[i];
        for(int j=1; j<4; j++){
            key.push_back(j);
            register int face_ind = l_cell_to_face[key];
            vector<int> key1;
            key1.push_back(face_ind);
            double deltx;
            double delty;
            double nx; 
            double ny;
            if(is_boundary_face(face_ind)){           
                vector<double> face_coor = face_center(face_ind);
                deltx = face_coor[0] - center_coor[0];
                delty = face_coor[1] - center_coor[1];
                key1.push_back(1);
                nx = sn[key1];
                key1.erase(key1.end()-1);
                key1.push_back(2);
                ny = sn[key1];
                if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
                    nx*=(-1);
                }
                if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
                    ny*=(-1);
                }
                if((deltx != 0) and (delty != 0)){
                    bc += (area(face_ind))*P_faces[face_ind-1]*((nx/deltx) +(ny/delty));
                    ac += (area(face_ind))*((nx/deltx) + (ny/delty));
                    gDiff.push_back((-1)*(area(face_ind))*((nx/deltx) + (ny/delty)));
                }
                else if (deltx == 0){
                    bc += (area(face_ind))*P_faces[face_ind-1]*(ny/delty);
                    ac += (area(face_ind))*(ny/delty);
                    gDiff.push_back((-1)*(area(face_ind))*(ny/delty));
                }else if (delty == 0){
                    bc += (area(face_ind))*P_faces[face_ind-1]*(nx/deltx);
                    ac += (area(face_ind))*(nx/deltx);
                    gDiff.push_back((-1)*(area(face_ind))*(nx/deltx));
                }
                key1.clear();
            }else{
                key1.push_back(1);
                nx = sn[key1];
                register int out_center_ind = l_face_to_cell[key1]; 
                key1.erase(key1.end()-1);
                key1.push_back(2);
                vector<double> c2 = centers[out_center_ind];
                ny = sn[key1];
                if(out_center_ind == i){
                    out_center_ind = l_face_to_cell[key1];
                    c2 = centers[out_center_ind];
                }
                deltx = c2[0] - center_coor[0];
                delty = c2[1] - center_coor[1];
                if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
                    nx*=(-1);
                }
                if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
                    ny*=(-1);
                }
                if((deltx != 0) and (delty != 0)){
                    gauss_seidel_coeff[i-1][out_center_ind-1] = (-1)*(area(face_ind))*((nx/deltx) +(ny/delty));
                    ac += (area(face_ind))*((nx/deltx) +(ny/delty));
                    gDiff.push_back((-1)*(area(face_ind))*((nx/deltx) +(ny/delty)));
                }else if (deltx == 0){
                    gauss_seidel_coeff[i-1][out_center_ind-1] = (-1)*(area(face_ind))*(ny/delty);
                    ac += (area(face_ind))*(ny/delty);
                    gDiff.push_back((-1)*(area(face_ind))*(ny/delty));
                }else if (delty ==0){
                    gauss_seidel_coeff[i-1][out_center_ind-1] = (-1)*(area(face_ind))*(nx/deltx);
                    ac += (area(face_ind))*(nx/deltx);
                    gDiff.push_back((-1)*(area(face_ind))*(nx/deltx));
                }
                key1.clear();
            }
            key.erase(key.end()-1);
        }
        gauss_seidel_coeff[i-1][i-1] = ac;
        gauss_seidel_b[i-1] = bc;
        ac = 0;
        bc = 0; 
        key.clear();
        gDiff.clear();
    }
    gauss_seidel();
}
void nodes_coordinates_generator(){
    vector<double> v1; 
    int cnt = 1;
    for(double i=0; i<n_side; i++){
        for(double j=0; j<n_side; j++){
            v1.push_back(j*step);
            v1.push_back(i*step);
            nodes[cnt] = v1;
            v1.clear();
            cnt++;
        }
    }
}
void cell_to_face(){
    vector<int> v; 
    v.push_back(n*n+1); 
    v.push_back(n*(2*n+1)+1); 
    v.push_back(1); 
    vector<int> cell_to_face_temp; 
    for(int i = 1; i < n * n + 1; i++){ 
        if (i % n == 1 and i > 1){ 
            v.at(1) += 1; 
        } 
        cell_to_face_temp.push_back(i); 
        for(int j = 1; j < 4; j++){ 
            cell_to_face_temp.push_back(j); 
            l_cell_to_face[cell_to_face_temp] = v[j-1]; 
            cell_to_face_temp.erase(cell_to_face_temp.end()-1); 
        } 
        cell_to_face_temp.clear(); 
        v.at(0) += 1; 
        v.at(1) += 1; 
        v.at(2) += 1; 
    } 
    v.at(0) = 1; 
    v.at(1) = n*(n+1)+1; 
    v.at(2) = n*(2*n+1)+2; 
    for(int i = n*n+1; i < 2*n*n+1; i++){ 
        if (i % n == 1 and i != n*n+1){ 
            v.at(2) += 1; 
        } 
        cell_to_face_temp.push_back(i); 
        for(int j = 1; j < 4; j++){ 
            cell_to_face_temp.push_back(j); 
            l_cell_to_face[cell_to_face_temp] = v[j-1]; 
            cell_to_face_temp.erase(cell_to_face_temp.end()-1); 
        } 
        cell_to_face_temp.clear(); 
        v.at(0) += 1; 
        v.at(1) += 1; 
        v.at(2) += 1; 
    } 
}
void face_to_cell(){
    vector<int> v; 
    v.push_back(1); 
    v.push_back(n*n+1); 
    vector<int> face_to_cell_temp; 
    for(int i = 1; i < n * n + 1; i++){ 
        face_to_cell_temp.push_back(i); 
        for(int j = 1; j < 3; j++){ 
            face_to_cell_temp.push_back(j); 
            l_face_to_cell[face_to_cell_temp] = v[j-1]; 
            face_to_cell_temp.erase(face_to_cell_temp.end()-1); 
        } 
        face_to_cell_temp.clear(); 
        v.at(0) += 1; 
        v.at(1) += 1; 
    } 
    v.at(0) = 1; 
    v.at(1) = 1; 
    for(int i = n * n + 1; i < n*(n+1)+1; i++){ 
        face_to_cell_temp.push_back(i); 
        for(int j = 1; j < 3; j++){ 
            face_to_cell_temp.push_back(j); 
            l_face_to_cell[face_to_cell_temp] = v[j-1]; 
            face_to_cell_temp.erase(face_to_cell_temp.end()-1); 
        } 
        face_to_cell_temp.clear(); 
        v.at(0) += 1; 
        v.at(1) += 1; 
    } 
    v.at(0) = n*n+1; 
    v.at(1) = n+1; 
    for(int i = n*(n+1)+1; i < 2*n*n+1; i++){ 
        face_to_cell_temp.push_back(i); 
        for(int j = 1; j < 3; j++){ 
            face_to_cell_temp.push_back(j); 
            l_face_to_cell[face_to_cell_temp] = v[j-1]; 
            face_to_cell_temp.erase(face_to_cell_temp.end()-1); 
        } 
        face_to_cell_temp.clear(); 
        v.at(0) += 1; 
        v.at(1) += 1; 
    } 
    v.at(0) = 2*n*n-n+1; 
    v.at(1) = 2*n*n-n+1; 
    for(int i = 2*n*n+1; i < n*(2*n+1)+1; i++){ 
        face_to_cell_temp.push_back(i); 
        for(int j = 1; j < 3; j++){ 
            face_to_cell_temp.push_back(j); 
            l_face_to_cell[face_to_cell_temp] = v[j-1]; 
            face_to_cell_temp.erase(face_to_cell_temp.end()-1); 
        } 
        face_to_cell_temp.clear(); 
        v.at(0) += 1; 
        v.at(1) += 1; 
    } 
    v.at(0) = n*n; 
    v.at(1) = 1; 
    for(int i = n*(2*n+1)+1; i < n*(3*n+2)+1; i++){ 
        face_to_cell_temp.push_back(i); 
        if(i%(n+1) == 2){ 
            for(int j = 1; j < 3; j++){ 
                face_to_cell_temp.push_back(j); 
                l_face_to_cell[face_to_cell_temp] = v[1]; 
                face_to_cell_temp.erase(face_to_cell_temp.end()-1);  
            } 
            v.at(0) += 1; 
            v.at(1) += 1; 
        }else if(i%(n+1) == 1){ 
            for(int j = 1; j < 3; j++){ 
                face_to_cell_temp.push_back(j); 
                l_face_to_cell[face_to_cell_temp] = v[0]; 
                face_to_cell_temp.erase(face_to_cell_temp.end()-1); 
            } 
        } 
        else{ 
            for(int j = 1; j < 3; j++){ 
                face_to_cell_temp.push_back(j); 
                l_face_to_cell[face_to_cell_temp] = v[j-1]; 
                face_to_cell_temp.erase(face_to_cell_temp.end()-1); 
                 
            } 
            v.at(0) += 1; 
            v.at(1) += 1; 
        } 
        face_to_cell_temp.clear(); 
    }
}
void face_to_node(){
    vector<int> v;
    v.push_back(n+2); 
    v.push_back(2); 
    vector<int> face_to_node_temp; 
    for(int i = 1; i < n * n + 1; i++){ 
        if (i % n == 1 and i > 1){ 
            v.at(0) += 1; 
            v.at(1) += 1; 
        } 
        face_to_node_temp.push_back(i); 
        for(int j = 1; j < 3; j++){ 
            face_to_node_temp.push_back(j); 
            l_face_to_node[face_to_node_temp] = v[j-1]; 
            face_to_node_temp.erase(face_to_node_temp.end()-1); 
        } 
        face_to_node_temp.clear(); 
        v.at(0) += 1; 
        v.at(1) += 1; 
    } 
    v.at(0) = 2; 
    v.at(1) = 1; 
    for(int i = n*n+1; i < n*(2*n+1)+1; i++){ 
        if (i % n == 1 and i != n*n+1){ 
            v.at(0) += 1; 
            v.at(1) += 1; 
        } 
        face_to_node_temp.push_back(i); 
        for(int j = 1; j < 3; j++){ 
            face_to_node_temp.push_back(j); 
            l_face_to_node[face_to_node_temp] = v[j-1]; 
            face_to_node_temp.erase(face_to_node_temp.end()-1); 
        } 
        face_to_node_temp.clear(); 
        v.at(0) += 1; 
        v.at(1) += 1; 
    } 
    v.at(0) = 1; 
    v.at(1) = n+2; 
    for(int i =  n*(2*n+1)+1; i < n*(3*n+2)+1; i++){ 
        face_to_node_temp.push_back(i); 
        for(int j = 1; j < 3; j++){ 
            face_to_node_temp.push_back(j); 
            l_face_to_node[face_to_node_temp] = v[j-1]; 
            face_to_node_temp.erase(face_to_node_temp.end()-1); 
        } 
        face_to_node_temp.clear(); 
        v.at(0) += 1; 
        v.at(1) += 1; 
    }  
}
void l_node_to_cell(){
    vector<int> v;
    v.push_back(2); 
    v.push_back(n*n+1); 
    v.push_back(n+1); 
    v.push_back(n*(n+1)+1); 
    v.push_back(n+2); 
    v.push_back(n*n+2); 
    vector<int> node_to_cell_temp; 
    for(int i = n+3; i < n*(n+1); i++){ 
        if (i % (n+1) == 1 or i % (n+1) == 0){ 
            continue; 
        }else{ 
            node_to_cell_temp.push_back(i); 
            if (i % (n+1) == 2 and i != n+3){ 
                v.at(0) += 1; 
                v.at(1) += 1; 
                v.at(2) += 1; 
                v.at(3) += 1; 
                v.at(4) += 1; 
                v.at(5) += 1; 
            } 
            for(int j = 1; j < 7; j++){ 
                node_to_cell_temp.push_back(j); 
                node_to_cell[node_to_cell_temp] = v[j-1]; 
                node_to_cell_temp.erase(node_to_cell_temp.end()-1);  
            } 
            node_to_cell_temp.clear(); 
            v.at(0) += 1; 
            v.at(1) += 1; 
            v.at(2) += 1; 
            v.at(3) += 1; 
            v.at(4) += 1; 
            v.at(5) += 1; 
        } 
    }
}
void center_coordinates_generator(){
    vector<int> v;
    vector<double> v1;
    vector<int> v1_;
    vector<vector<double> > points_;
    for(int i = 1; i < center_num + 1; i++){
        v.push_back(i);
        v.push_back(1);
        register int face1 = l_cell_to_face[v];
        v.erase(v.end()-1);
        v.push_back(2);
        register int face2 = l_cell_to_face[v];
        v.erase(v.end()-1);
        v.push_back(3);
        register int face3 = l_cell_to_face[v];
        v1_.push_back(face1);
        v1_.push_back(1);
        vector<double> node1 = nodes[l_face_to_node[v1_]];
        points_.push_back(node1);
        v1_.erase(v1_.end()-1);
        v1_.push_back(2);
        vector<double> node2 = nodes[l_face_to_node[v1_]];
        points_.push_back(node2);
        v1_.clear();
        v1_.push_back(face2);
        v1_.push_back(1);
        vector<double> node3 = nodes[l_face_to_node[v1_]];
        if(count(points_.begin(), points_.end(), node3) == 0){
            points_.push_back(node3);
        }
        v1_.erase(v1_.end()-1);
        v1_.push_back(2);
        vector<double> node4 = nodes[l_face_to_node[v1_]];
        if(count(points_.begin(), points_.end(), node4) == 0){
            points_.push_back(node4);
        }
        v1.push_back((points_[0][0] + points_[1][0] + points_[2][0]) / 3);
        v1.push_back((points_[0][1] + points_[1][1] + points_[2][1]) / 3);
        centers[i]=v1;
        points_.clear();
        v1_.clear();
        v.clear();
        v1.clear();
    }
}

void set_boundary_conditions(){
    register int bound1 = 1;//bottom
    register int bound2 = n_side;
    register int bound[n_side];
    for(int i=0; i<n_side; i++){
        bound[i] = bound1+i;
    }
    register int bound3[n];
    for(int i=0; i<n; i++){
        bound3[i] = i+1+center_num/2;
    }//bottom
    // int bound1 = n_side;//left
    // int bound2 = node_num; 
    // int bound[n_side];
    // for(int i=0; i<n_side; i++){
    //     bound[i] = (i+1)*n_side;
    // }
    // int bound3[n];
    // bound3[0] = (center_num/2 + n_side*n)+1+n; 
    // for(int i=1; i<n; i++){
    //     bound3[i] = bound3[i-1] + n_side;
    // }//left 
    // int bound1 = 1;//right
    // int bound2 = node_num - n; 
    // int bound[n_side];
    // bound[0] =1;
    // for(int i=1; i<n_side; i++){
    //     bound[i] = (i*n_side+1);
    // }
    // int bound3[n];
    // bound3[0] = (center_num/2 + n_side*n)+1; 
    // for(int i=1; i<n; i++){
    //     bound3[i] = bound3[i-1] + n_side;
    // }//right 
    // int bound1 = node_num - n; //top 
    // int bound2 = node_num; 
    // int bound[n_side];
    // for(int i=0; i<n_side; i++){
    //     bound[i] = bound1+i;
    // }
    // int bound3[n];
    // for(int i=0; i<n; i++){
    //     bound3[i] = (center_num/2 + n_side*n)-i;
    // }//top
    for(int i=0; i<node_num; i++){
        P_nodes[i] = 0;
    } 
    for(int i=0; i<n_side; i++){
        P_nodes[bound[i]-1] = 1;
    } 
    for(int i=0; i<face_num; i++){
        P_faces[i]=0;
    }
    for(int i=0; i<n; i++){
        P_faces[bound3[i]-1] = 1;
    } 
    for(int i=0; i<center_num; i++){
        P_centers[i]=0;
        P_centers_copy[i]=0;
    } 
}
int main(){
    unsigned int start_time =  clock();
    nodes_coordinates_generator();
    cell_to_face();
    face_to_cell();
    face_to_node();
    l_node_to_cell();
    center_coordinates_generator();
    vector<double> v1;
    vector<int> v1_;
    for(int i=1; i<face_num+1; i++){
        v1 = normal(i);
        v1_.push_back(i);
        v1_.push_back(1);
        sn[v1_] = v1[0];
        v1_.erase(v1_.end()-1);
        v1_.push_back(2);
        sn[v1_] = v1[1];
        v1_.clear();
        v1.clear();
    }
    set_boundary_conditions();
    fvm();
    for(int i=0; i<center_num; i++){
        cout<<P_centers[i]<<" ";
        if((i+1)%n == 0){
            cout<<endl;
        }
    }
    for(int i = n_side+1; i < node_num-n_side; i++){
        if (i % n_side != 1 and i % n_side != 0){
            P_nodes[i-1] = cell_to_vertex_interpolation(i);
        }
    }
    for(int i=0; i<node_num; i++){
        cout<<P_nodes[i]<<" ";
        if((i+1)%n_side == 0){
            cout<<endl;
        }
    }
    unsigned int end_time = clock(); 
    unsigned int search_time = end_time - start_time;
    cout<<"overall: "<<search_time/CLOCKS_PER_SEC<<"sec"<<endl;
    return 0;
}