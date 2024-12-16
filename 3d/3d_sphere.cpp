//FVM in 3d
#include <bits/stdc++.h>  
using namespace std;
const int node_num=718;
map<int, vector<double> > nodes;
const int center_num=3242;
map<int, vector<double> > centers;
int face_num;
map<vector<int> , int> l_cell_to_face;
map<vector<int> , int> l_cell_to_node;
map<vector<int> , int> l_face_to_cell;
map<vector<int> , int> l_face_to_node;
map<int, vector<int> > node_to_cell;
map<vector<int> , double> sn;
vector<int> boundary_face;
vector<int> nodes_boundary;
double P_nodes[node_num];
vector<double> P_faces;
double P_centers[center_num];
double P_centers_copy[center_num];
const double e=0.000001;
double gauss_seidel_coeff[center_num][center_num];
double gauss_seidel_b[center_num];
vector<double> split_double(string line, char delimeter){//checked
    string my_str = "";
    vector<double> v;
    for(int i; i<line.size(); i++){
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
vector<int> split_int(string line, char delimeter){//checked
    string my_str = "";
    vector<int> v;
    for(int i; i<line.size(); i++){
        if((line[i] == delimeter) or (i==line.size()-1)){
            if(line[i+1] == delimeter){
                continue;
            }
            if(i==line.size()-1){
                my_str+=line[i];  
            }
            if(my_str[0]== 'C'){
                my_str = "";
                continue;
            }
            char* temp = new char[my_str.length()+1];
            strcpy(temp, my_str.c_str());
            int to_insert = atof(temp);
            v.push_back(to_insert);
            my_str = "";
        }else{
            my_str+=line[i];  
        }
    }
    return v;
}
void nodes_coor(){//checked
    string myText;
    ifstream NodesFile("sphere_coor.txt");
    vector<double> v; 
    // int cnt = 1;
    // int cnt_node = 1;
    while (getline(NodesFile, myText)) {
        // if(cnt%3 ==0){
        v = split_double(myText, ' ');
        int cnt = v[0]; 
        v.erase(v.begin());
        nodes[cnt] = v;
        // cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<endl;
        v.clear();
        // cnt_node++;
        // }
        // cnt++;
        // if(cnt==node_num){
        //     break;
        // }
    }
    NodesFile.close();
}
void cell_to_node(){//Checked
    string myText;
    ifstream NodesFile("sphere_links.txt");
    vector<int> v; 
    vector<int> v1; 
    while (getline(NodesFile, myText)) {
        v = split_int(myText, ' ');
        int cell = v[0];
        v.erase(v.begin());
        v.erase(v.begin());
        for(int i=1; i<5; i++){
            v1.push_back(cell);
            v1.push_back(i);
            l_cell_to_node[v1] = v[i-1];
            if(node_to_cell.find(v[i-1]) != node_to_cell.end()){
               node_to_cell[v[i-1]].push_back(cell);
            }else{
                vector<int> temp;
                temp.push_back(cell);
                node_to_cell[v[i-1]] = temp;
                temp.clear();
            }
            v1.clear();
        }
        v.clear();
    }
    NodesFile.close();
}
void center_coordinates_generator(){//Checked without output
    vector<int> v;
    vector<double> v1;
    for(int i = 1; i < center_num + 1; i++){
        v.push_back(i);
        v.push_back(1);
        register int node1 = l_cell_to_node[v];
        v.erase(v.end()-1);
        v.push_back(2);
        register int node2 = l_cell_to_node[v];
        v.erase(v.end()-1);
        v.push_back(3);
        register int node3 = l_cell_to_node[v];
        v.erase(v.end()-1);
        v.push_back(4);
        register int node4 = l_cell_to_node[v];
        v1.push_back((nodes[node1][0] + nodes[node2][0] + nodes[node3][0]+ nodes[node4][0]) / 4);
        v1.push_back((nodes[node1][1] + nodes[node2][1] + nodes[node3][1]+ nodes[node4][1]) / 4);
        v1.push_back((nodes[node1][2] + nodes[node2][2] + nodes[node3][2]+ nodes[node4][2]) / 4);
        centers[i]=v1;
        v.clear();
        v1.clear();
    }
}
void face_links(){//boundary_faces are not ordered
    map< set<int>, vector<int> > nodes_to_cell;
    vector<int> key;
    vector<int> key1;
    register int cnt = 1;
    vector<int> trg_nodes;
    set<int> nodes_ind;
    vector< vector<int> > c_t_f(center_num);
    for(int i=1; i<center_num+1; i++){
        key.push_back(i);
        for(int j=1; j<5; j++){
            key.push_back(j);
            register int trg_node = l_cell_to_node[key];
            trg_nodes.push_back(trg_node);
            key.erase(key.end()-1);
        }
        key.clear();
        int arr[4][3] = {{0,1,2}, {0,1,3}, {0,2,3}, {1,2,3}};
        for(int j=1; j<5; j++){
            nodes_ind.insert(trg_nodes[arr[j-1][0]]);
            nodes_ind.insert(trg_nodes[arr[j-1][1]]); 
            nodes_ind.insert(trg_nodes[arr[j-1][2]]);
            if(nodes_to_cell.find(nodes_ind) == nodes_to_cell.end()){
                key1.push_back(i);
                nodes_to_cell[nodes_ind] = key1; 
                key1.clear(); 
            }else{
                nodes_to_cell[nodes_ind].push_back(i);
            }
            nodes_ind.clear();           
        }
        trg_nodes.clear();
    }
    face_num = nodes_to_cell.size();
    map< set<int>, vector<int> >::iterator m_it=nodes_to_cell.begin();
    register int f_cnt = 1;
    while(m_it!=nodes_to_cell.end()){
        key1 = m_it->second; 
        set <int>::iterator s_it = m_it->first.begin();
        register int s_cnt = 1;
        while(s_it != m_it->first.end()){
            key.push_back(f_cnt); 
            key.push_back(s_cnt);
            l_face_to_node[key] = *s_it;
            s_cnt++;
            s_it++;
            key.clear();
        }
        key.push_back(f_cnt);
        for(int i=1; i<key1.size()+1; i++){
            c_t_f[key1[i-1]-1].push_back(f_cnt);
            key.push_back(i);
            l_face_to_cell[key] = key1[i-1]; 
            key.erase(key.end()-1);
            if(key1.size() == 1){
                boundary_face.push_back(f_cnt);
                key.push_back(i+1);
                key.erase(key.end()-1);
                l_face_to_cell[key] = key1[i-1]; 
            }
            
        }
        f_cnt++;
        m_it++;
        key1.clear();
        key.clear();
    }
    // key.push_back(boundary_face[1]);
    // key.push_back(1);
    // key1.push_back(boundary_face[2]);
    // key1.push_back(0);
    // if(l_face_to_node[key] != l_face_to_node[key1]){
    //     int temp = boundary_face[0];
    //     boundary_face[0] = boundary_face[1];
    //     boundary_face[1] = temp;
    // }
    key.clear();
    key1.clear();
    for(int i=1; i<center_num+1; i++){
        key.push_back(i);
        for(int j=1; j<c_t_f[i-1].size()+1; j++){
            key.push_back(j);
            l_cell_to_face[key] = c_t_f[i-1][j-1];
            key.erase(key.end()-1);
        }
        key.clear();
    }
}
vector<double> face_center(int face_ind){//Check changes in num node
    vector <int> key; 
    key.push_back(face_ind); 
    key.push_back(1);
    vector<double> n1 = nodes[l_face_to_node[key]];
    key.erase(key.end()-1);
    key.push_back(2);
    vector<double> n2 = nodes[l_face_to_node[key]];
    key.erase(key.end()-1);
    key.push_back(3);
    vector<double> n3 = nodes[l_face_to_node[key]];
    vector<double> ans; 
    ans.push_back((n1[0] + n2[0] + n3[0])/3);
    ans.push_back((n1[1] + n2[1] + n3[1])/3);
    ans.push_back((n1[2] + n2[2] + n3[2])/3);
    return ans;
}
double find_max(double  v[center_num]){
    double max = v[1];
    for(int i=1; i<center_num; i++){
        if(v[i]>max){
            max = v[i];
        }
    }
    return max;
}
void gauss_seidel(){
    double to_check = 1;
    while(abs(to_check - find_max(P_centers))>e){
        to_check = find_max(P_centers);
        for(int i=0; i<center_num; i++){
            P_centers_copy[i] = gauss_seidel_b[i];
            for(int j=0; j<center_num; j++){
                if (j == i){
                    continue;
                }
                P_centers_copy[i] = P_centers_copy[i] - (gauss_seidel_coeff[i][j] * P_centers[j]);
                P_centers[i] = (P_centers_copy[i]/gauss_seidel_coeff[i][i]);
                // cout<<gauss_seidel_coeff[i][i]<<endl;
                // cout << P_centers[j] << endl;
            }
        }
    }
}
double min_ = 0; 
vector<double> normal(int face_ind){
    vector <double> ans; 
    vector <int> key; 
    key.push_back(face_ind); 
    key.push_back(1);
    // cout<<"x: "<<l_face_to_node[key]<<endl;

    vector<double> n1 = nodes[l_face_to_node[key]];
    key.erase(key.end()-1);
    key.push_back(2);
    // cout<<"y: "<<l_face_to_node[key]<<endl;
    vector<double> n2 = nodes[l_face_to_node[key]];
    key.erase(key.end()-1); 
    key.push_back(3);
    // cout<<"z: "<<l_face_to_node[key]<<endl;
    vector<double> n3 = nodes[l_face_to_node[key]];
    // if(face_ind == 1855){
        // cout<<"x "<<n1[0]<<" "<<n2[0]<<" "<<n3[0]<<endl;
        // cout<<"y "<<n1[1]<<" "<<n2[1]<<" "<<n3[1]<<endl;
        // cout<<"z "<<n1[2]<<" "<<n2[2]<<" "<<n3[2]<<endl;
    // }
 

    vector<double> v1;
    vector<double> v2;
    for(int i=0; i<3; i++){
        v1.push_back(n2[i] - n1[i]);
        v2.push_back(n2[i] - n3[i]);
    }
    // cout<<"v1 "<<v1[0]<<" "<<v1[1]<<" "<<v1[2]<<endl;
    // cout<<"v2 "<<v2[0]<<" "<<v2[1]<<" "<<v2[2]<<endl;

    ans.push_back((v1[1]*v2[2])-(v1[2]*v2[1]));
    ans.push_back((v1[2]*v2[0])-(v1[0]*v2[2]));
    ans.push_back((v1[0]*v2[1])-(v2[0]*v1[1]));
    // if(ans[0] == ans[1] == ans[2] == 0){
    //     cout<<ans[0]<<" "<<ans[1]<<" "<<ans[2]<<endl;
    // }
    if(ans[0] <min_){
        min_ = ans[0];
    }
    if(ans[1] <min_){
        min_ = ans[1];
    }
    if(ans[2] <min_){
        min_ = ans[2];
    }
    
    double length=0;
    for(int i=0; i<3; i++){
        length += pow(ans[i], 2);
        // if(length == 0){
        //     cout<<length<<endl;
        // }
    }
    // if(length == 0){
    //     cout<<length<<endl;
    //     cout<<ans[0]<<" "<<ans[1]<<" "<<ans[2]<<endl;
    //     cout<<face_ind<<endl;

    // }
    length = (pow(length, 0.5));
    for(int i=0; i<3; i++){
        ans[i] = ans[i]/length;
    }
    return ans; 
}
double area(int face_ind){
    vector <double> ans; 
    vector <int> key; 
    key.push_back(face_ind); 
    key.push_back(1);
    vector<double> n1 = nodes[l_face_to_node[key]];
    // cout << "n1 " << n1[0] << " " << n1[1] << " " << n1[2] << endl;
    key.erase(key.end()-1);
    key.push_back(2);
    vector<double> n2 = nodes[l_face_to_node[key]];
    // cout << "n2 " << n2[0] << " " << n2[1] << " " << n2[2] << endl;
    key.erase(key.end()-1); 
    key.push_back(3);
    vector<double> n3 = nodes[l_face_to_node[key]];
    // cout << "n3 " << n3[0] << " " << n3[1] << " " << n3[2] << endl;
    vector<double> v1;
    vector<double> v2;
    for(int i=0; i<3; i++){
        v1.push_back(n2[i] - n1[i]);
        v2.push_back(n2[i] - n3[i]);
    }
    // cout << "v1 " << v1[0] << " " << v1[1] << " " << v1[2] << endl;
    // cout << "v2 " << v2[0] << " " << v2[1] << " " << v2[2] << endl;
    ans.push_back((v1[1]*v2[2])-(v1[2]*v2[1]));
    ans.push_back((v1[2]*v2[0])-(v1[0]*v2[2]));
    ans.push_back((v1[0]*v2[1])-(v2[0]*v1[1]));
    
    double length=0;
    for(int i=0; i<3; i++){
        length += pow(ans[i], 2);
    }
    // cout<<"length "<<length<<endl;
    length = (pow(length, 0.5));
    return length/2; 
}
double cell_to_vertex_interpolation(int node_ind){ //check
    vector<double> node_coor = nodes[node_ind]; 
    vector<double> center_coor;  
    int cell_ind; 
    double nom = 0; 
    double denom = 0; 
    for(int i=0; i<node_to_cell[node_ind].size(); i++){ 
        cell_ind = node_to_cell[node_ind][i]; 
        // cout << node_ind << " " << i << " " << cell_ind << endl;
        center_coor = centers[cell_ind]; 
        double d = pow((pow((center_coor[0] - node_coor[0]), 2)   +   pow((center_coor[1] - node_coor[1]), 2) + pow((center_coor[2] - node_coor[2]), 2)), 0.5); 
        double P = P_centers[cell_ind - 1]; 
        // cout << "d = " << d << endl;
        // cout << "P = " << P << endl;
        nom += P/d;
         
        denom += 1/d;  
        // cout << nom << " " << denom << endl;
    } 
    
    return nom/denom; 
}
void fvm(){//check
    for(int i=1; i<center_num+1; i++){
        vector<double> gDiff;
        vector<int> key;
        key.push_back(i);
        double ac=0;
        double bc=0;
        vector<double> center_coor = centers[i];
        for(int j=1; j<5; j++){
            key.push_back(j);
            register int face_ind = l_cell_to_face[key];
            vector<int> key1;
            key1.push_back(face_ind);
            double deltx;
            double delty;
            double deltz;
            double nx; 
            double ny;
            double nz;
            double valid_denom = 0;
            if(find(boundary_face.begin(), boundary_face.end(), face_ind) != boundary_face.end()){   
                vector<double> face_coor = face_center(face_ind);
                deltx = face_coor[0] - center_coor[0];
                delty = face_coor[1] - center_coor[1];
                deltz = face_coor[2] - center_coor[2];
                key1.push_back(1);
                nx = sn[key1];
                key1.erase(key1.end()-1);
                key1.push_back(2);
                ny = sn[key1];
                key1.erase(key1.end()-1);
                key1.push_back(3);
                nz = sn[key1];
                if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
                    nx*=(-1);
                }
                if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
                    ny*=(-1);
                }
                if(((deltz>0) and (nz<0)) or ((deltz<0) and (nz>0))){
                    nz*=(-1);
                }
                valid_denom = 0;
                if(deltx != 0){
                    valid_denom += (nx/deltx);
                }
                if(delty != 0){
                    valid_denom += (ny/delty);
                }
                if(deltz != 0){
                    valid_denom += (nz/deltz);
                }
                // if(isnan(valid_denom)){
                //     cout<<"deltx: "<<deltx<<endl;
                //     cout<<"delty: "<<delty<<endl;
                //     cout<<"deltz: "<<deltz<<endl;
                //     cout<<"nx: "<<nx<<endl;
                //     cout<<"ny: "<<ny<<endl;
                //     cout<<"nz: "<<nz<<endl;
                // }
                bc += (area(face_ind))*P_faces[face_ind-1]*valid_denom;
                ac += (area(face_ind))*valid_denom;
                gDiff.push_back((-1)*(area(face_ind))*valid_denom);
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
                key1.erase(key1.end()-1);
                key1.push_back(3);
                nz = sn[key1];
                deltx = c2[0] - center_coor[0];
                delty = c2[1] - center_coor[1];
                deltz = c2[2] - center_coor[2];
                if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
                    nx*=(-1);
                }
                if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
                    ny*=(-1);
                }
                if(((deltz>0) and (nz<0)) or ((deltz<0) and (nz>0))){
                    nz*=(-1);
                }
                valid_denom = 0;
                if(deltx != 0){
                    valid_denom += (nx/deltx);
                }
                if(delty != 0){
                    valid_denom += (ny/delty);
                }
                if(deltz != 0){
                    valid_denom += (nz/deltz);
                }
                // if(isnan(valid_denom)){
                //     cout<<"deltx: "<<deltx<<endl;
                //     cout<<"delty: "<<delty<<endl;
                //     cout<<"deltz: "<<deltz<<endl;
                //     cout<<"nx: "<<nx<<endl;
                //     cout<<"ny: "<<ny<<endl;
                //     cout<<"nz: "<<nz<<endl;
                // }
                gauss_seidel_coeff[i-1][out_center_ind-1] = (-1)*(area(face_ind))*valid_denom;
                // cout<<"a: "<<gauss_seidel_coeff[i-1][out_center_ind-1]<<endl;
                ac += (area(face_ind))*valid_denom;
                // cout<< "AREA " << area(face_ind) << endl; 
                gDiff.push_back((-1)*(area(face_ind))*valid_denom);
                key1.clear();
            }
            key.erase(key.end()-1);
        }
        gauss_seidel_coeff[i-1][i-1] = ac;
        // cout<<"ac: "<<gauss_seidel_coeff[i-1][i-1]<<endl;
        gauss_seidel_b[i-1] = bc;
        // cout<<"bc: "<<gauss_seidel_b[i-1]<<endl;
        ac = 0;
        bc = 0; 
        key.clear();
        gDiff.clear();
    }
    gauss_seidel();
}
void set_boundary_conditions(){//rewrite
    register int bound1 = 100;
    register int bound2 = 340;
    for(int i=0; i<face_num; i++){
        P_faces.push_back(0);
    }
    for(int i=0; i<node_num; i++){
        P_nodes[i] = 0;
    }
    vector<int> v; 
    for(int i=bound1; i<bound2; i++){
        P_faces[boundary_face[i]-1] = 1;
        v.push_back(boundary_face[i]); 
        v.push_back(1);
        P_nodes[l_face_to_node[v]-1] = 1;
        // cout<<nodes[l_face_to_node[v]][0]<<" "<<nodes[l_face_to_node[v]][1]<<" "<<nodes[l_face_to_node[v]][2]<<endl;
        nodes_boundary.push_back(l_face_to_node[v]);
        v.erase(v.end()-1);
        v.push_back(2);
        P_nodes[l_face_to_node[v]-1] = 1;
        // cout<<nodes[l_face_to_node[v]][0]<<" "<<nodes[l_face_to_node[v]][1]<<" "<<nodes[l_face_to_node[v]][2]<<endl;
        nodes_boundary.push_back(l_face_to_node[v]);
        v.erase(v.end()-1);
        v.push_back(3);
        P_nodes[l_face_to_node[v]-1] = 1; 
        // cout<<nodes[l_face_to_node[v]][0]<<" "<<nodes[l_face_to_node[v]][1]<<" "<<nodes[l_face_to_node[v]][2]<<endl;
        nodes_boundary.push_back(l_face_to_node[v]);
        v.clear();
    }
    // bound1 = 300;
    // bound2 = 450;
    // for(int i=bound1; i<bound2; i++){
    //     P_faces[boundary_face[i]-1] = 1;
    //     v.push_back(boundary_face[i]); 
    //     v.push_back(1);
    //     P_nodes[l_face_to_node[v]-1] = 1;
    //     // cout<<nodes[l_face_to_node[v]][0]<<" "<<nodes[l_face_to_node[v]][1]<<" "<<nodes[l_face_to_node[v]][2]<<endl;
    //     nodes_boundary.push_back(l_face_to_node[v]);
    //     v.erase(v.end()-1);
    //     v.push_back(2);
    //     P_nodes[l_face_to_node[v]-1] = 1; 
    //     // cout<<nodes[l_face_to_node[v]][0]<<" "<<nodes[l_face_to_node[v]][1]<<" "<<nodes[l_face_to_node[v]][2]<<endl;
    //     nodes_boundary.push_back(l_face_to_node[v]);
    //     v.erase(v.end()-1);
    //     v.push_back(3);
    //     P_nodes[l_face_to_node[v]-1] = 1;
    //     // cout<<nodes[l_face_to_node[v]][0]<<" "<<nodes[l_face_to_node[v]][1]<<" "<<nodes[l_face_to_node[v]][2]<<endl;
    //     nodes_boundary.push_back(l_face_to_node[v]);
    //     v.clear();
    // }
    for(int i=0; i<center_num; i++){
        P_centers[i]=0;
        P_centers_copy[i]=0;
    } 
}
int main(){
    nodes_coor();
    cout<<"1"<<endl;
    cout<<nodes.size()<<endl;
    cell_to_node();
    cout<<"2"<<endl;

    center_coordinates_generator();
    cout<<"3"<<endl;

    face_links();
    cout<<"4"<<endl;

    vector<double> v1;
    vector<int> v1_;
    cout<<face_num<<endl;
    for(int i=1; i<face_num+1; i++){
        v1 = normal(i);
        // cout<<v1[0]<<" "<<v1[1]<<" "<<v1[2]<<endl;
        if(isnan(v1[0])){
            cout<<"face "<<i<<endl;
            cout<<"nx "<<v1[0]<<endl;
        }
        if(isnan(v1[1])){
            cout<<"face "<<i<<endl;
            cout<<"ny "<<v1[1]<<endl;
        }
        if(isnan(v1[2])){
            cout<<"face "<<i<<endl;
            cout<<"nz "<<v1[2]<<endl;
        }
        v1_.push_back(i);
        v1_.push_back(1);
        sn[v1_] = v1[0];
        v1_.erase(v1_.end()-1);
        v1_.push_back(2);
        sn[v1_] = v1[1];
        v1_.erase(v1_.end()-1);
        v1_.push_back(3);
        sn[v1_] = v1[2];
        v1_.clear();
        v1.clear();
    }
    cout<<sn.size()<<endl;
    cout<<"min: "<<min_<<endl;
    set_boundary_conditions();
    cout<<"5"<<endl;

    // cout<<boundary_face.size()<<endl;
    fvm();
    cout<<"6"<<endl;

    for(int i=1; i<node_num+1; i++){
        if(find(nodes_boundary.begin(), nodes_boundary.end(), i) != nodes_boundary.end()){
                continue;
        }else{
            P_nodes[i-1] = cell_to_vertex_interpolation(i);
        }
    }
    cout<<boundary_face.size()<<endl;
    ofstream MyFile("pressure_3d_sphere.dat");
    MyFile << "VARIABLES= \"X\", \"Y\", \"Z\",\"P\"" <<endl;
    for(int i=0; i<node_num; i++){
        MyFile << nodes[i+1][0] << " " << nodes[i+1][1] << " " << nodes[i+1][2] <<" "<<P_nodes[i]<<endl;
    }
    // MyFile.close();
    // for(int i=0; i<center_num; i++){
    //     MyFile << centers[i+1][0] << " " << centers[i+1][1] << " " << centers[i+1][2] <<" "<<P_centers[i]<<endl;
    // }
    MyFile.close();
    sort(P_centers, P_centers+center_num, greater<double>());
    ofstream MFile("pressure_sorted_sphere.dat");
    // MFile << "VARIABLES= \"X\", \"Y\", \"Z\",\"P\"" <<endl;
    for(int i=0; i<center_num; i++){
        MFile <<P_centers[i]<<endl;
    }
    MFile.close();
    return 0;
}   