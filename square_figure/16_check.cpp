// Triangle mesh generator 
//NODES numeration -1 
#include <bits/stdc++.h>  
using namespace std;

const int n = 9;
const int n_side = n+1;
const int node_num=n_side*n_side;
const double step = 4.5/double(n);
map<int, vector<double> > nodes;
const int center_num=(n_side-1)*(n_side-1)*2;
map<int, vector<double> > centers;
const int face_num = ((n_side-1)*(n_side)*2) + (center_num/2);

map<vector<int> , int> l_cell_to_face;
map<vector<int> , int> l_face_to_cell;
map<vector<int> , int> l_face_to_node;
map<vector<int> , int> node_to_cell;
map<vector<int> , double> sn;


map<int, double> P_nodes;
map<int, double>P_faces;
map<int, double> P_centers;
map<int, double> P_centers_copy;

const double e=0.000001;

double gauss_seidel_coeff[center_num][center_num];
double gauss_seidel_b[center_num];

double find_max(map<int, double>  v){
    double max = v[1];
    for(int i=1; i<v.size()+1; i++){
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
    //cout << "node: " << node_ind << endl;
    //cout << "**************************" << endl;
    //cout << "node_ind = " << node_ind << endl; 
    vector<double> node_coor = nodes[node_ind]; 
    vector<int> key; 
    key.push_back(node_ind); 
    key.push_back(1); 
    int cell1 = node_to_cell[key]; 
    key.erase(key.end()-1);  
    key.push_back(2); 
    int cell2 = node_to_cell[key]; 
    key.erase(key.end()-1);  
    key.push_back(3); 
    int cell3 = node_to_cell[key]; 
    key.erase(key.end()-1);  
    key.push_back(4); 
    int cell4 = node_to_cell[key]; 
    key.erase(key.end()-1);  
    key.push_back(5); 
    int cell5 = node_to_cell[key];
    key.erase(key.end()-1);  
    key.push_back(6); 
    int cell6 = node_to_cell[key];  
    //cout << "Cells" << " " << cell1 << " " << cell2 << " " << cell3 << " " << cell4 << " " << cell5 << " " << cell6 << endl;
    vector<double> center1 = centers[cell1]; 
    vector<double> center2 = centers[cell2]; 
    vector<double> center3 = centers[cell3]; 
    vector<double> center4 = centers[cell4]; 
    vector<double> center5 = centers[cell5]; 
    vector<double> center6 = centers[cell6];
    double P1 = P_centers[cell1]; 
    double P2 = P_centers[cell2]; 
    double P3 = P_centers[cell3]; 
    double P4 = P_centers[cell4]; 
    double P5 = P_centers[cell5]; 
    double P6 = P_centers[cell6];
    double d1 = pow(pow((center1[0] - node_coor[0]), 2) + pow((center1[1] - node_coor[1]), 2), 0.5); 
    double d2 = pow(pow((center2[0] - node_coor[0]), 2) + pow((center2[1] - node_coor[1]), 2), 0.5); 
    double d3 = pow(pow((center3[0] - node_coor[0]), 2) + pow((center3[1] - node_coor[1]), 2), 0.5); 
    double d4 = pow(pow((center4[0] - node_coor[0]), 2) + pow((center4[1] - node_coor[1]), 2), 0.5); 
    double d5 = pow(pow((center5[0] - node_coor[0]), 2) + pow((center5[1] - node_coor[1]), 2), 0.5); 
    double d6 = pow(pow((center6[0] - node_coor[0]), 2) + pow((center6[1] - node_coor[1]), 2), 0.5); 
    // cout << "center1: " <<  center1[0] << " " << center1[1] << endl;
    // cout << "center2: " <<  center2[0] << " " << center2[1] << endl;
    // cout << "center3: " <<  center3[0] << " " << center3[1] << endl;
    // cout << "center4: " <<  center4[0] << " " << center4[1] << endl;
    // cout << "center5: " <<  center5[0] << " " << center5[1] << endl;
    // cout << "center6: " <<  center6[0] << " " << center6[1] << endl;
    // cout << "P " << P1 << " " << P2 << " " << P3 << " " << P4 << " " << P5 << " " << P6 << endl;
    // cout << "d " << d1 << " " << d2 << " " << d3 << " " << d4 << " " << d5 << " " << d6 << endl;
    return (P1/d1 + P2/d2 + P3/d3 + P4/d4 + P5/d5 + P6/d6)/(1/d1 + 1/d2 + 1/d3 + 1/d4 + 1/d5 + 1/d6); 
}

void gauss_seidel(){

    double to_check = 1;
 
    while(abs(to_check - find_max(P_centers))>e){
        to_check = find_max(P_centers);
        for(int i=1; i<center_num+1; i++){
        
            P_centers_copy[i] = (gauss_seidel_b[i-1] / gauss_seidel_coeff[i-1][i-1]);

            for(int j=1; j<center_num+1; j++){
                if (j == i){
                    continue;
                }
                P_centers_copy[i] = P_centers_copy[i] - ((gauss_seidel_coeff[i-1][j-1] / gauss_seidel_coeff[i-1][i-1]) * P_centers[j]);

                P_centers[i] = P_centers_copy[i];
            }
        }
    }
}

void fvm(int i){
    // for(int i=1; i<center_num+1; i++){
        vector<double> gDiff;
        vector<int> key;
        key.push_back(i);
        double ac=0;
        double bc=0;
        vector<double> center_coor = centers[i];
        for(int j=1; j<4; j++){
            key.push_back(j);
            int face_ind = l_cell_to_face[key];
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
                    bc += (area(face_ind))*P_faces[face_ind]*((nx/deltx) +(ny/delty));
                    ac += (area(face_ind))*((nx/deltx) + (ny/delty));
                    gDiff.push_back((-1)*(area(face_ind))*((nx/deltx) + (ny/delty)));
                }
                else if (deltx == 0){
                    bc += (area(face_ind))*P_faces[face_ind]*(ny/delty);
                    ac += (area(face_ind))*(ny/delty);
                    gDiff.push_back((-1)*(area(face_ind))*(ny/delty));
                }else if (delty == 0){
                    bc += (area(face_ind))*P_faces[face_ind]*(nx/deltx);
                    ac += (area(face_ind))*(nx/deltx);
                    gDiff.push_back((-1)*(area(face_ind))*(nx/deltx));
                }
                    
                cout<<"face = "<<face_ind<<endl;
                cout<<"deltx = "<<deltx<<endl;
                cout<<"delty = "<<delty<<endl;
                cout<<"nx = "<<nx<<endl;
                cout<<"ny = "<<ny<<endl;
                cout<<"area = "<<area(face_ind)<<endl;
                cout<<"P = "<<P_faces[face_ind]<<endl;
                cout<<"   ac = "<<ac<<endl;
                cout<<"   a = "<<gDiff[gDiff.size()-1]<<endl;
                cout<<"   bc = "<<bc<<endl;

                key1.clear();

            }else{
               
                key1.push_back(1);
                nx = sn[key1];
                int out_center_ind = l_face_to_cell[key1]; 
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
               
                cout<<"face = "<<face_ind<<endl;
                cout<<"outcenter = "<<out_center_ind<<endl;
                cout<<"deltx = "<<deltx<<endl;
                cout<<"delty = "<<delty<<endl;
                cout<<"nx = "<<nx<<endl;
                cout<<"ny = "<<ny<<endl;
                cout<<"area = "<<area(face_ind)<<endl;
                cout<<"P = "<<P_faces[face_ind]<<endl;
                cout<<"   ac = "<<ac<<endl;
                cout<<"   a = "<<gDiff[gDiff.size()-1]<<endl;
                cout<<"   bc = "<<bc<<endl;
                key1.clear();
            }
            key.erase(key.end()-1);
        }
        cout<<"ac = "<<ac<<endl;
        cout<<"bc = "<<bc<<endl;
        for(int m=1; m<gDiff.size()+1; m++){
            cout<<"a"<<m<<" = "<<gDiff[m-1]<<endl;
        }
        gauss_seidel_coeff[i-1][i-1] = ac;
        gauss_seidel_b[i-1] = bc;
        ac = 0;
        bc = 0; 
        key.clear();
        gDiff.clear();
        for(int s=0; s<center_num; s++){
            cout<<gauss_seidel_coeff[i][s]<<" ";
        }
        cout<<endl;
    // }
    // for(int s=0; s<center_num; s++){
    // for(int l=0; l<center_num; l++){

    //         cout<<gauss_seidel_coeff[s][l]<<" ";
    //     }
    // }
    // gauss_seidel();
    // jacobi();
}

int main(){
    // cout<<face_num<<" :face num"<<endl;
    // cout<<center_num<<" :center num"<<endl;
    // cout<<node_num<<" :node num"<<endl;


    //Nodes coordinates generator
    vector<double> v; 
    int cnt = 1;
    for(double i=0; i<n_side; i++){
        for(double j=0; j<n_side; j++){
            v.push_back(j*step);
            v.push_back(i*step);
            nodes[cnt] = v;
            v.clear();
            cnt++;
        }
    }

    //Nodes Outputing 

    map<int, vector<double> >::iterator it_node=nodes.begin(); 
    while(it_node!=nodes.end()){ 
        // pair<vector<int> , int> res=*it; 
        // cout<<res.first.first<<" "<<res.first.second<<" = "<<res.second<<endl;// the same --> cout<<(*it)first<<"--> "<<(*it)second<<endl; 
        cout<<it_node->first<<"-->"<<it_node->second[0]<<" "<<it_node->second[1]<<endl; 
        it_node++; 
    } 


    //Cell to face generator
    vector<int> l; 
    l.push_back(n*n+1); 
    l.push_back(n*(2*n+1)+1); 
    l.push_back(1); 
    vector<int> cell_to_face_temp; 
    // cout << "============" << endl; 
    for(int i = 1; i < n * n + 1; i++){ 
        if (i % n == 1 and i > 1){ 
            l.at(1) += 1; 
        } 
        cell_to_face_temp.push_back(i); 
        for(int j = 1; j < 4; j++){ 
            cell_to_face_temp.push_back(j); 
            l_cell_to_face[cell_to_face_temp] = l[j-1]; 
            cell_to_face_temp.erase(cell_to_face_temp.end()-1); 
        } 
        cell_to_face_temp.clear(); 
        l.at(0) += 1; 
        l.at(1) += 1; 
        l.at(2) += 1; 
    } 
    l.at(0) = 1; 
    l.at(1) = n*(n+1)+1; 
    l.at(2) = n*(2*n+1)+2; 
    for(int i = n*n+1; i < 2*n*n+1; i++){ 
        if (i % n == 1 and i != n*n+1){ 
            l.at(2) += 1; 
        } 
        cell_to_face_temp.push_back(i); 
        for(int j = 1; j < 4; j++){ 
            cell_to_face_temp.push_back(j); 
            l_cell_to_face[cell_to_face_temp] = l[j-1]; 
            cell_to_face_temp.erase(cell_to_face_temp.end()-1); 
        } 
        cell_to_face_temp.clear(); 
        l.at(0) += 1; 
        l.at(1) += 1; 
        l.at(2) += 1; 
    } 
     
     
    map<vector<int> , int>::iterator it=l_cell_to_face.begin(); 
    while(it!=l_cell_to_face.end()){ 
        // pair<vector<int> , int> res=*it; 
        // cout<<res.first.first<<" "<<res.first.second<<" = "<<res.second<<endl;// the same --> cout<<(*it)first<<"--> "<<(*it)second<<endl; 
        cout<<it->first[0]<<" "<<it->first[1]<<"-->"<<it->second<<endl; 
        it++; 
    } 


    //Face to cell generator
    vector<int> v1; 
    for(int i=1; i<n_side+1; i++){
        v1.push_back(i);
        v1.push_back(1);
        l_face_to_cell[v1]=i;
        v1.erase(v1.end()-1);
        v1.push_back(2);
        l_face_to_cell[v1]=i;
        v1.clear();
    }
    int stopper = (n_side*2) -1;
    int start = n_side-1;
    int cell1 = n_side-1;
    int cell2 = 1;
    int temp;
    for(int i=n_side+1; i<face_num-n+1; i++){
        if(i == (start + stopper)){
            if(stopper == ((n_side*2) -1)){
                start = i+1;
                temp = cell1;
                cell1=cell2; 
                cell2 = temp+1;
                stopper = n_side-1;
                v1.push_back(i);
                v1.push_back(1);
                l_face_to_cell[v1]=cell1;
                v1.erase(v1.end()-1);
                v1.push_back(2);
                l_face_to_cell[v1]=cell1;
                v1.clear();
                cell1 = cell1+1;
                cell2= cell1-(n_side-1);
            }else if (stopper == n_side-1){
                start = i-1;
                stopper = ((n_side*2) -1);
                v1.push_back(i);
                v1.push_back(1);
                l_face_to_cell[v1]=cell2;
                v1.erase(v1.end()-1);
                v1.push_back(2);
                l_face_to_cell[v1]=cell2;
                v1.clear();
                cell2 = cell2; 
                cell1 = cell1 -1;  
            }
            continue;
        }

        if(stopper == (n_side*2) -1){
            temp = cell1;
            cell1=cell2; 
            cell2 = temp+1;
            v1.push_back(i);
            v1.push_back(1);
            l_face_to_cell[v1]=cell1;
            v1.erase(v1.end()-1);
            v1.push_back(2);
            l_face_to_cell[v1]=cell2;
            v1.clear();
        }else if(stopper == n_side-1){
            v1.push_back(i);
            v1.push_back(1);
            l_face_to_cell[v1]=cell1;
            v1.erase(v1.end()-1);
            v1.push_back(2);
            l_face_to_cell[v1]=cell2;
            v1.clear();
            cell1++;
            cell2++;
        }

    }
    for(int i=face_num-n+1; i<face_num+1; i++){
        v1.push_back(i);
        v1.push_back(1);
        l_face_to_cell[v1]= cell2;
        v1.erase(v1.end()-1);
        v1.push_back(2);
        l_face_to_cell[v1]= cell2;
        cell2++;
        v1.clear();
    }
    // map<vector<int> , int>::iterator it=l_face_to_cell.begin();
    // while(it!=l_face_to_cell.end()){
    //     cout<<it->first[0]<<" "<<it->first[1]<<"-->"<<it->second<<endl;
    //     it++;
    // }

    //Face to node generator
    vector<int> l1; 
    l1.push_back(n+2); 
    l1.push_back(2); 
    vector<int> face_to_node_temp; 
    for(int i = 1; i < n * n + 1; i++){ 
        if (i % n == 1 and i > 1){ 
            l1.at(0) += 1; 
            l1.at(1) += 1; 
        } 
        face_to_node_temp.push_back(i); 
        for(int j = 1; j < 3; j++){ 
            face_to_node_temp.push_back(j); 
            l_face_to_node[face_to_node_temp] = l1[j-1]; 
            face_to_node_temp.erase(face_to_node_temp.end()-1); 
        } 
        face_to_node_temp.clear(); 
        l1.at(0) += 1; 
        l1.at(1) += 1; 
    } 
     
    l1.at(0) = 2; 
    l1.at(1) = 1; 
    for(int i = n*n+1; i < n*(2*n+1)+1; i++){ 
        if (i % n == 1 and i != n*n+1){ 
            l1.at(0) += 1; 
            l1.at(1) += 1; 
        } 
        face_to_node_temp.push_back(i); 
        for(int j = 1; j < 3; j++){ 
            face_to_node_temp.push_back(j); 
            l_face_to_node[face_to_node_temp] = l1[j-1]; 
            face_to_node_temp.erase(face_to_node_temp.end()-1); 
        } 
        face_to_node_temp.clear(); 
        l1.at(0) += 1; 
        l1.at(1) += 1; 
    } 
     
    l1.at(0) = 1; 
    l1.at(1) = n+2; 
    for(int i =  n*(2*n+1)+1; i < n*(3*n+2)+1; i++){ 
        face_to_node_temp.push_back(i); 
        for(int j = 1; j < 3; j++){ 
            face_to_node_temp.push_back(j); 
            l_face_to_node[face_to_node_temp] = l1[j-1]; 
            face_to_node_temp.erase(face_to_node_temp.end()-1); 
        } 
        face_to_node_temp.clear(); 
        l1.at(0) += 1; 
        l1.at(1) += 1; 
    } 
     
     
     
    // map<vector<int> , int>::iterator it1=l_face_to_node.begin(); 
    // while(it1!=l_face_to_node.end()){ 
    //     // pair<vector<int> , int> res=*it; 
    //     // cout<<res.first.first<<" "<<res.first.second<<" = "<<res.second<<endl;// the same --> cout<<(*it)first<<"--> "<<(*it)second<<endl; 
    //     cout<<it1->first[0]<<" "<<it1->first[1]<<"-->"<<it1->second<<endl; 
    //     it1++; 
    // } 
    // cout << "============" << endl;

    //Node to cell generator 
    //cout << "Node to cell" << endl;
    vector<int> l2; 
    l2.push_back(2); 
    l2.push_back(n*n+1); 
    l2.push_back(n+1); 
    l2.push_back(n*(n+1)+1); 
    l2.push_back(n+2); 
    l2.push_back(n*n+2); 
    vector<int> node_to_cell_temp; 
    for(int i = n+3; i < n*(n+1); i++){ 
        if (i % (n+1) == 1 or i % (n+1) == 0){ 
            continue; 
        }else{ 
            node_to_cell_temp.push_back(i); 
            if (i % (n+1) == 2 and i != n+3){ 
                l2.at(0) += 1; 
                l2.at(1) += 1; 
                l2.at(2) += 1; 
                l2.at(3) += 1; 
                l2.at(4) += 1; 
                l2.at(5) += 1; 
            } 
            for(int j = 1; j < 7; j++){ 
                node_to_cell_temp.push_back(j); 
                node_to_cell[node_to_cell_temp] = l2[j-1]; 
                node_to_cell_temp.erase(node_to_cell_temp.end()-1); 
                 
            } 
            node_to_cell_temp.clear(); 
            l2.at(0) += 1; 
            l2.at(1) += 1; 
            l2.at(2) += 1; 
            l2.at(3) += 1; 
            l2.at(4) += 1; 
            l2.at(5) += 1; 
        } 
    }
    // map<vector<int> , int>::iterator it2=node_to_cell.begin(); 
    // while(it2!=node_to_cell.end()){ 
    //     // pair<vector<int> , int> res=*it; 
    //     // cout<<res.first.first<<" "<<res.first.second<<" = "<<res.second<<endl;// the same --> cout<<(*it)first<<"--> "<<(*it)second<<endl; 
    //     cout<<it2->first[0]<<" "<<it2->first[1]<<"-->"<<it2->second<<endl; 
    //     it2++; 
    // } 



    //Center Coordinate generator 
    vector<int> v_;
    vector<int> v1_;
    vector<vector<double> > points_;
    // cout << "-------------------" << endl;
    
    for(int i = 1; i < center_num + 1; i++){
        v_.push_back(i);
        v_.push_back(1);
        int face1 = l_cell_to_face[v_];
        v_.erase(v_.end()-1);
        v_.push_back(2);
        int face2 = l_cell_to_face[v_];
        v_.erase(v_.end()-1);
        v_.push_back(3);
        int face3 = l_cell_to_face[v_];
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
        v.push_back((points_[0][0] + points_[1][0] + points_[2][0]) / 3);
        v.push_back((points_[0][1] + points_[1][1] + points_[2][1]) / 3);
        centers[i]=v;
        points_.clear();
        v1_.clear();
        v_.clear();
        v.clear();
    }
    //Centers Outputing 
    // map<int, vector<double> >::iterator it_center=centers.begin(); 
    // while(it_center!=centers.end()){ 
    //     cout<<it_center->first<<"-->"<<it_center->second[0]<<" "<<it_center->second[1]<<endl; 
    //     it_center++; 
    // } 

    //Surface normal filling

    for(int i=1; i<face_num+1; i++){
        v = normal(i);
        v1.push_back(i);
        v1.push_back(1);
        sn[v1] = v[0];
        v1.erase(v1.end()-1);
        v1.push_back(2);
        sn[v1] = v[1];
        v1.clear();
        v.clear();
    }

    //Pressure filling
    
    // int bound1 = 1;//bottom
    // int bound2 = n_side;
    // int bound[n_side];
    // for(int i=0; i<n_side; i++){
    //     bound[i] = bound1+i;
    // }
    // int bound3[n];
    // for(int i=0; i<n; i++){
    //     bound3[i] = i+1;
    // }//bottom

    // int bound1 = n_side;//left
    // int bound2 = node_num; 
    // int bound[n_side];
    // for(int i=0; i<n_side; i++){
    //     bound[i] = (i+1)*n_side;
    // }
    // int bound3[n];
    // for(int i=0; i<n; i++){
    //     bound3[i] = (n_side + (n*2))*(i+1);
    // }//left 

    // int bound1 = 1;//right
    // int bound2 = node_num - n; 
    // int bound[n_side];
    // bound[0] =1;
    // for(int i=1; i<n_side; i++){
    //     bound[i] = (i*n_side+1);
    // }
    // int bound3[n];
    // bound3[0] = n_side; 
    // for(int i=1; i<n; i++){
    //     bound3[i] = (n_side + (n*2)) + bound3[i-1];
    // }//right 

    int bound1 = node_num - n; //top 
    int bound2 = node_num; 
    int bound[n_side];
    for(int i=0; i<n_side; i++){
        bound[i] = bound1+i;
    }
    int bound3[n];
    for(int i=0; i<n; i++){
        bound3[i] = face_num - i;
    }//top

    // cout<<"bound1: "<<bound1<<endl;

    // for(int i=0; i<n_side; i++){
    //     cout<<bound[i]<<endl;
    // }

    // cout<<"bound2: "<<bound2<<endl;
    // for(int i=0; i<n; i++){
    //     cout<<bound3[i]<<endl;
    // }

        //Nodes Filling
    for(int i=0; i<node_num; i++){
        P_nodes[i] = 0;
    } 

    for(int i=0; i<n_side; i++){
        P_nodes[bound[i]] = 1;
    } 
        //Faces Filling
    for(int i=0; i<face_num; i++){
        P_faces[i]=0;
    }

    for(int i=0; i<n; i++){
        P_faces[bound3[i]] = 1;
    } 

        //Centers filling
    for(int i=1; i<center_num+1; i++){
        P_centers[i]=0;
        P_centers_copy[i]=0;

    } 
    // fvm(1);
    //Pressure in centers outputting
    // for(int i=1; i<center_num+1; i++){
    //     cout<<P_centers[i]<<" ";
    //     if(i%n == 0){
    //         cout<<endl;
    //     }
    // }

    // //Cell to vertex intepolation
    // for(int i = n_side+1; i < node_num-n_side; i++){
    //     if (i % n_side != 1 and i % n_side != 0){
    //         P_nodes[i] = cell_to_vertex_interpolation(i);
    //     //cout << P_nodes[i] << endl;
    //     }
    // }

    // //Pressure in nodes outputting
    // for(int i=1; i<node_num+1; i++){
    //     cout<<P_nodes[i]<<" ";
    //     if(i%n_side == 0){
    //         cout<<endl;
    //     }
    // }

    return 0;
}